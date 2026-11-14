//! Tau-leaping of the voxel RDME (opt-in, `ssa.tau_leap = true`).
//!
//! The channels of a leap are the reactions of every voxel and the hops of
//! every species through every face (`X(v) -> X(w)`, a first-order channel
//! of propensity `n D / h^2`). As in the unit cell ([`crate::ssa::tau`]):
//!
//! - A channel is critical when the population it consumes is below
//!   [`RdmeTau::n_critical`] (for a reaction: fewer firings left than that).
//!   Critical channels fire one at a time.
//! - The leap bounds the relative change of every propensity by `epsilon`:
//!   the bound of Cao, Gillespie and Petzold is taken over every (voxel,
//!   species) population, with the hops in and out in its mean and variance.
//! - Every non-critical channel fires a Poisson number of times. A leap that
//!   would make a count negative is rejected and redone with half the step.
//! - A leap reads every voxel, so it only pays when it replaces many events:
//!   when fewer than `10 + 4 x voxels` non-critical events are expected, the
//!   run takes exact next-subvolume events instead, in batches that double
//!   while leaping keeps being refused.
//! - A leap never crosses the time it is asked to advance to (seed and
//!   snapshot times), and the coarsening test runs on the leaped events too.

use crate::ssa::rdme::{Level, Rdme};
use crate::ssa::tau::{Order, species_orders};
use crate::voxel::grid::VoxelGrid;

/// Scratch and counters of the voxel tau-leaping.
#[derive(Debug, Clone)]
pub struct RdmeTau {
    /// Bound on the relative change of a propensity over a leap.
    pub epsilon: f64,
    /// A channel whose consumed population is below this is critical.
    pub n_critical: i64,
    /// Leaps taken.
    pub n_leaps: u64,
    /// Leaps rejected for a negative count.
    pub n_rejected: u64,
    /// Exact events fired between leaps.
    pub n_exact: u64,
    /// Size of the last batch of exact events taken when leaping was refused.
    exact_batch: u64,
    /// Exact events still to fire before leaping is considered again. It is
    /// kept across calls: a sparse grid then pays for one leap decision per
    /// batch, not one per output time.
    exact_left: u64,
    orders: Vec<Order>,
    mu: Vec<f64>,
    sigma2: Vec<f64>,
    /// Total critical rate of every voxel.
    critical: Vec<f64>,
    trial: Vec<i64>,
}

impl RdmeTau {
    fn new(rdme: &Rdme<'_>, epsilon: f64) -> Self {
        Self {
            epsilon,
            n_critical: 10,
            n_leaps: 0,
            n_rejected: 0,
            n_exact: 0,
            exact_batch: 0,
            exact_left: 0,
            orders: species_orders(rdme.setup.model()),
            mu: Vec::new(),
            sigma2: Vec::new(),
            critical: Vec::new(),
            trial: Vec::new(),
        }
    }
}

/// Rates of one leap decision.
struct Plan {
    /// Total rate of every channel.
    a0: f64,
    /// Total rate of the critical channels.
    a_critical: f64,
    /// The leap bounding the change of every propensity.
    tau: f64,
}

/// Call `f(axis, neighbour)` for every face of voxel `v`.
fn for_each_face(grid: &VoxelGrid, v: usize, mut f: impl FnMut(usize, usize)) {
    let c = grid.coords(v);
    let strides = [1, grid.dims[0], grid.dims[0] * grid.dims[1]];
    for a in 0..3 {
        if c[a] > 0 {
            f(a, v - strides[a]);
        }
        if c[a] + 1 < grid.dims[a] {
            f(a, v + strides[a]);
        }
    }
}

/// Whether reaction `j` (positive propensity) has fewer than `n_critical`
/// firings left for the counts `n`.
fn reaction_is_critical(lv: &Level, j: usize, n: &[i64], n_critical: i64) -> bool {
    lv.model.reactions[j]
        .need
        .iter()
        .any(|(s, c)| n[*s] / c < n_critical)
}

/// A critical channel chosen to fire.
enum Fire {
    Reaction { v: usize, j: usize },
    Hop { v: usize, w: usize, s: usize },
}

impl Rdme<'_> {
    /// Tau-leaping with the accuracy `epsilon` (`None`: exact).
    pub fn set_tau_leaping(&mut self, epsilon: Option<f64>) {
        self.tau = epsilon.map(|e| Box::new(RdmeTau::new(self, e)));
    }

    /// Counters of the tau-leaping (leaps, rejected leaps, exact events).
    pub fn tau_stats(&self) -> Option<(u64, u64, u64)> {
        self.tau
            .as_ref()
            .map(|t| (t.n_leaps, t.n_rejected, t.n_exact))
    }

    /// Rates, critical channels and the leap bound for the current counts.
    fn plan_leap(&self, tau: &mut RdmeTau) -> Plan {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let (n_sp, n_vox) = (setup.n_sp, lv.grid.n_voxels());
        tau.mu.clear();
        tau.mu.resize(n_vox * n_sp, 0.0);
        tau.sigma2.clear();
        tau.sigma2.resize(n_vox * n_sp, 0.0);
        tau.critical.clear();
        tau.critical.resize(n_vox, 0.0);
        let (mut a0, mut a_critical) = (0.0, 0.0);
        for v in 0..n_vox {
            let n = &self.counts[v * n_sp..(v + 1) * n_sp];
            let mut critical = 0.0;
            for (j, r) in lv.model.reactions.iter().enumerate() {
                let a = lv.model.propensity(j, n);
                if a <= 0.0 {
                    continue;
                }
                a0 += a;
                if reaction_is_critical(lv, j, n, tau.n_critical) {
                    critical += a;
                    continue;
                }
                for &(s, d) in &r.delta {
                    let d = d as f64;
                    tau.mu[v * n_sp + s] += d * a;
                    tau.sigma2[v * n_sp + s] += d * d * a;
                }
            }
            for s in 0..n_sp {
                let (count, d) = (n[s], setup.diffusion[s]);
                if count <= 0 || d <= 0.0 {
                    continue;
                }
                let is_critical = count < tau.n_critical;
                for_each_face(&lv.grid, v, |axis, w| {
                    let a = count as f64 * d * lv.inv_h2[axis];
                    a0 += a;
                    if is_critical {
                        critical += a;
                    } else {
                        tau.mu[v * n_sp + s] -= a;
                        tau.sigma2[v * n_sp + s] += a;
                        tau.mu[w * n_sp + s] += a;
                        tau.sigma2[w * n_sp + s] += a;
                    }
                });
            }
            tau.critical[v] = critical;
            a_critical += critical;
        }
        let mut bound = f64::INFINITY;
        for v in 0..n_vox {
            for s in 0..n_sp {
                let order = tau.orders[s];
                if order == Order::None && setup.diffusion[s] <= 0.0 {
                    continue;
                }
                let i = v * n_sp + s;
                let (mu, sigma2) = (tau.mu[i], tau.sigma2[i]);
                if mu == 0.0 && sigma2 <= 0.0 {
                    continue;
                }
                let x = self.counts[i];
                let allowed = (tau.epsilon * x as f64 / order.g(x)).max(1.0);
                if mu != 0.0 {
                    bound = bound.min(allowed / mu.abs());
                }
                if sigma2 > 0.0 {
                    bound = bound.min(allowed * allowed / sigma2);
                }
            }
        }
        Plan {
            a0,
            a_critical,
            tau: bound,
        }
    }

    /// The critical channel at `target` in `[0, total critical rate)`.
    fn pick_critical(&self, tau: &RdmeTau, target: f64) -> Option<Fire> {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let n_sp = setup.n_sp;
        let mut acc = 0.0;
        let mut voxel = None;
        for (v, c) in tau.critical.iter().enumerate() {
            if *c > 0.0 {
                voxel = Some((v, acc));
                acc += c;
                if target < acc {
                    break;
                }
            }
        }
        let (v, before) = voxel?;
        let n = &self.counts[v * n_sp..(v + 1) * n_sp];
        let mut acc = before;
        let mut last = None;
        for j in 0..lv.model.reactions.len() {
            let a = lv.model.propensity(j, n);
            if a > 0.0 && reaction_is_critical(lv, j, n, tau.n_critical) {
                acc += a;
                last = Some(Fire::Reaction { v, j });
                if target < acc {
                    return last;
                }
            }
        }
        for s in 0..n_sp {
            let (count, d) = (n[s], setup.diffusion[s]);
            if count <= 0 || d <= 0.0 || count >= tau.n_critical {
                continue;
            }
            let mut found = None;
            for_each_face(&lv.grid, v, |axis, w| {
                if found.is_some() {
                    return;
                }
                acc += count as f64 * d * lv.inv_h2[axis];
                last = Some(Fire::Hop { v, w, s });
                if target < acc {
                    found = Some(Fire::Hop { v, w, s });
                }
            });
            if found.is_some() {
                return found;
            }
        }
        // Round-off left the target past the last partial sum.
        last
    }

    /// Draw the firings of a step `dt` into `tau.trial`. Returns the
    /// reactions and hops fired, or `None` if a count would go negative.
    fn draw_leap(&mut self, tau: &mut RdmeTau, dt: f64, fire: Option<Fire>) -> Option<(u64, u64)> {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let (n_sp, n_vox) = (setup.n_sp, lv.grid.n_voxels());
        tau.trial.clear();
        tau.trial.extend_from_slice(&self.counts);
        let (mut reactions, mut hops) = (0u64, 0u64);
        for v in 0..n_vox {
            let n = &self.counts[v * n_sp..(v + 1) * n_sp];
            for (j, r) in lv.model.reactions.iter().enumerate() {
                let a = lv.model.propensity(j, n);
                if a <= 0.0 || reaction_is_critical(lv, j, n, tau.n_critical) {
                    continue;
                }
                let k = self.rng.poisson(a * dt);
                if k > 0 {
                    reactions += k as u64;
                    for &(s, d) in &r.delta {
                        tau.trial[v * n_sp + s] += d * k;
                    }
                }
            }
            for s in 0..n_sp {
                let (count, d) = (n[s], setup.diffusion[s]);
                if count < tau.n_critical || d <= 0.0 {
                    continue;
                }
                let rng = &mut self.rng;
                let trial = &mut tau.trial;
                for_each_face(&lv.grid, v, |axis, w| {
                    let k = rng.poisson(count as f64 * d * lv.inv_h2[axis] * dt);
                    if k > 0 {
                        hops += k as u64;
                        trial[v * n_sp + s] -= k;
                        trial[w * n_sp + s] += k;
                    }
                });
            }
        }
        match fire {
            Some(Fire::Reaction { v, j }) => {
                reactions += 1;
                for &(s, d) in &lv.model.reactions[j].delta {
                    tau.trial[v * n_sp + s] += d;
                }
            }
            Some(Fire::Hop { v, w, s }) => {
                hops += 1;
                tau.trial[v * n_sp + s] -= 1;
                tau.trial[w * n_sp + s] += 1;
            }
            None => {}
        }
        tau.trial
            .iter()
            .all(|n| *n >= 0)
            .then_some((reactions, hops))
    }

    /// Advance to `t_stop` with tau-leaping where it pays.
    pub(crate) fn advance_tau(&mut self, t_stop: f64) {
        let Some(mut tau) = self.tau.take() else {
            return;
        };
        while self.t < t_stop {
            if tau.exact_left > 0 {
                if self.queue_stale {
                    self.rebuild();
                }
                while tau.exact_left > 0 && self.step(t_stop) {
                    tau.exact_left -= 1;
                    tau.n_exact += 1;
                    self.count_events(1);
                }
                if tau.exact_left > 0 {
                    // No event left before `t_stop`.
                    self.t = t_stop;
                }
                continue;
            }
            let n_vox = self.grid().n_voxels() as u64;
            let plan = self.plan_leap(&mut tau);
            if plan.a0 <= 0.0 {
                self.t = t_stop;
                break;
            }
            // Non-critical events one leap would fire.
            let remaining = t_stop - self.t;
            let mut dt = plan.tau.min(remaining);
            if plan.a_critical > 0.0 {
                dt = dt.min(1.0 / plan.a_critical);
            }
            let expected = (plan.a0 - plan.a_critical) * dt;
            if expected < (10 + 4 * n_vox) as f64 {
                // Exact events; the batch doubles while leaping is refused.
                tau.exact_batch = (tau.exact_batch * 2).clamp(8 * n_vox, 512 * n_vox);
                tau.exact_left = tau.exact_batch;
                continue;
            }
            tau.exact_batch = 0;
            let mut tau_leap = plan.tau;
            loop {
                let remaining = t_stop - self.t;
                let tau_critical = if plan.a_critical > 0.0 {
                    self.rng.exp1() / plan.a_critical
                } else {
                    f64::INFINITY
                };
                let (dt, fire, to_stop) = if remaining <= tau_leap && remaining <= tau_critical {
                    (remaining, false, true)
                } else if tau_critical < tau_leap {
                    (tau_critical, true, false)
                } else {
                    (tau_leap, false, false)
                };
                let chosen = if fire {
                    let target = self.rng.uniform() * plan.a_critical;
                    self.pick_critical(&tau, target)
                } else {
                    None
                };
                let Some((reactions, hops)) = self.draw_leap(&mut tau, dt, chosen) else {
                    tau.n_rejected += 1;
                    tau_leap *= 0.5;
                    continue;
                };
                std::mem::swap(&mut self.counts, &mut tau.trial);
                self.t = if to_stop { t_stop } else { self.t + dt };
                self.n_reactions += reactions;
                self.n_hops += hops;
                self.queue_stale = true;
                tau.n_leaps += 1;
                self.count_events(reactions + hops);
                break;
            }
        }
        self.tau = Some(tau);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::matrices::SimMatrices;
    use crate::ssa::rdme::RdmeSetup;
    use crate::ssa::rng::SsaRng;

    const CHEMISTRY: &str = r#"
        [[reactions.k_reaction]]
        reaction = "H_r -> H2"
        k_value = 3.0e3

        [[reactions.k_reaction]]
        reaction = "e_aq + OH_r -> OH_minus"
        k_value = 3.0e10
    "#;

    fn setup(dims: [usize; 3], h: f64, levels: u32) -> (SimEnv, RdmeSetup) {
        let env = SimEnv::from_config(&Config::from_toml_str(CHEMISTRY, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new(dims, [h; 3], [0.0; 3]).unwrap();
        let s = RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], levels, false).unwrap();
        (env, s)
    }

    #[test]
    fn leaps_conserve_molecules_and_stop_on_time() {
        let (env, s) = setup([4, 2, 2], 1e-6, 0);
        let idx = |l: &str| env.species.index_of(l).unwrap();
        let (h, h2, e, oh, ohm) = (
            idx("H_r"),
            idx("H2"),
            idx("e_aq"),
            idx("OH_r"),
            idx("OH_minus"),
        );
        for stream in 0..10 {
            let mut r = Rdme::new(&s, 0.0, SsaRng::new(7, stream));
            r.set_tau_leaping(Some(0.03));
            r.inject(&[(0, h, 40_000), (0, e, 30_000), (15, oh, 50_000)]);
            let mut t = 1e-6;
            while t < 2e-3 {
                r.advance_to(t);
                assert_eq!(r.t, t);
                assert!(r.counts.iter().all(|n| *n >= 0));
                let n = r.totals();
                assert_eq!(n[h] + n[h2], 40_000);
                assert_eq!(n[e] + n[ohm], 30_000);
                assert_eq!(n[oh] + n[ohm], 50_000);
                t *= 2.0;
            }
            let (leaps, _, exact) = r.tau_stats().unwrap();
            assert!(leaps > 0, "leaps were taken");
            // Far fewer steps than events.
            assert!(
                leaps + exact < r.n_events() / 5,
                "{leaps} leaps, {exact} exact"
            );
        }
    }

    #[test]
    fn sparse_grids_run_exact_events() {
        // A few molecules per voxel: every channel is critical, leaping is refused.
        let (env, s) = setup([4, 4, 4], 10e-9, 0);
        let e = env.species.index_of("e_aq").unwrap();
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(1, 0));
        r.set_tau_leaping(Some(0.03));
        r.inject(&[(0, e, 50)]);
        r.advance_to(1e-6);
        assert_eq!(r.t, 1e-6);
        assert_eq!(r.totals()[e], 50);
        let (leaps, rejected, exact) = r.tau_stats().unwrap();
        assert_eq!((leaps, rejected), (0, 0));
        assert_eq!(exact, r.n_events());
    }

    #[test]
    fn leaping_and_coarsening_work_together() {
        let (env, s) = setup([4, 4, 4], 1e-6, 2);
        let e = env.species.index_of("e_aq").unwrap();
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(2, 0));
        r.set_tau_leaping(Some(0.03));
        r.set_coarsening(Some(3.0));
        r.inject(&[(0, e, 200_000)]);
        // h^2 / D = 0.2 ms per voxel.
        r.advance_to(2e-2);
        assert_eq!(r.totals()[e], 200_000);
        assert_eq!(r.level(), 2, "history {:?}", r.level_history);
        assert!(r.tau_stats().unwrap().0 > 0);
    }
}
