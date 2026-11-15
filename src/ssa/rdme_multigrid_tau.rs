//! Tau-leaping of the multi-grid RDME (`ssa.tau_leap = true` with
//! `[ssa.species_levels]`), the counterpart of [`crate::ssa::rdme_tau`].
//!
//! The populations are the counts of every species on the voxels of the
//! level it lives on. The channels are the reactions of every site (a voxel
//! of the reaction's level) and the hops of every species through every
//! face of its level. The rules of the single-level leap hold:
//!
//! - a channel is critical when it has fewer than [`MultigridTau::n_critical`]
//!   firings left (for a reaction, counted on the summed counts of its
//!   site; for a hop, molecules in the voxel) and fires one at a time;
//! - the leap is the Cao-Gillespie-Petzold bound over every population;
//! - non-critical channels fire a Poisson number of times, and a leap that
//!   would make a count negative is rejected and halved;
//! - leaping is taken only when it replaces more than `10 + 4 x sites`
//!   events, otherwise exact events in doubling batches.
//!
//! A reaction whose species all live on its level changes the counts of
//! its voxel. A cross-level reaction fires one firing at a time, with the
//! rules of the exact path: each consumed molecule of a finer reactant is
//! taken from a contained voxel in proportion to the counts at the start of
//! the leap, and the products go where the finest consumed molecule was.
//! Its contribution to the mean and variance of the populations is spread
//! over the contained voxels in the same proportions.

use crate::ssa::rdme_multigrid::Multigrid;
use crate::ssa::tau::{Order, species_orders};
use crate::voxel::grid::VoxelGrid;

/// Scratch and counters of the multi-grid tau-leaping.
#[derive(Debug, Clone)]
pub struct MultigridTau {
    /// Bound on the relative change of a propensity over a leap.
    pub epsilon: f64,
    /// A channel with fewer firings left than this is critical.
    pub n_critical: i64,
    /// Leaps taken.
    pub n_leaps: u64,
    /// Leaps rejected for a negative count.
    pub n_rejected: u64,
    /// Exact events fired between leaps.
    pub n_exact: u64,
    exact_batch: u64,
    exact_left: u64,
    orders: Vec<Order>,
    /// First population of every species (its voxels on its current level).
    pop_off: Vec<usize>,
    mu: Vec<f64>,
    sigma2: Vec<f64>,
    /// Total critical rate of every site.
    critical: Vec<f64>,
    trial: Vec<i64>,
    /// Scratch: where the products of one firing may land, `(level, voxel, weight)`.
    spots: Vec<(usize, usize, f64)>,
}

impl MultigridTau {
    pub(crate) fn new(mg: &Multigrid<'_>, epsilon: f64) -> Self {
        Self {
            epsilon,
            n_critical: 10,
            n_leaps: 0,
            n_rejected: 0,
            n_exact: 0,
            exact_batch: 0,
            exact_left: 0,
            orders: species_orders(mg.setup.base.model()),
            pop_off: Vec::new(),
            mu: Vec::new(),
            sigma2: Vec::new(),
            critical: Vec::new(),
            trial: Vec::new(),
            spots: Vec::new(),
        }
    }
}

/// Rates of one leap decision.
struct Plan {
    a0: f64,
    a_critical: f64,
    tau: f64,
}

/// A critical channel chosen to fire.
enum Fire {
    Reaction { l: usize, v: usize, j: usize },
    Hop { s: usize, v: usize, w: usize },
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

impl Multigrid<'_> {
    /// Tau-leaping with the accuracy `epsilon` (`None`: exact).
    pub fn set_tau_leaping(&mut self, epsilon: Option<f64>) {
        self.tau = epsilon.map(|e| Box::new(MultigridTau::new(self, e)));
    }

    /// Counters of the tau-leaping (leaps, rejected leaps, exact events).
    pub fn tau_stats(&self) -> Option<(u64, u64, u64)> {
        self.tau
            .as_ref()
            .map(|t| (t.n_leaps, t.n_rejected, t.n_exact))
    }

    /// Voxels of level `to` inside voxel `v` of level `from >= to`.
    fn inside(&self, v: usize, from: usize, to: usize) -> impl Iterator<Item = usize> + use<'_> {
        let levels = &self.setup.base.levels;
        let c = levels[from].grid.coords(v);
        let n = 1usize << (from - to);
        let fine = &levels[to].grid;
        (0..n * n * n).map(move |k| {
            let (i, j, l) = (k % n, (k / n) % n, k / (n * n));
            fine.index(c[0] * n + i, c[1] * n + j, c[2] * n + l)
        })
    }

    /// Whether every species changed by reaction `j` lives on its level `l`.
    fn same_level(&self, l: usize, j: usize) -> bool {
        let lay = self.lay();
        self.setup.base.levels[l].model.reactions[j]
            .delta
            .iter()
            .all(|(s, _)| lay.level[*s] == l)
    }

    /// Whether reaction `j` has fewer than `n_critical` firings left on the
    /// counts loaded in `buf`.
    fn is_critical(&self, l: usize, j: usize, n_critical: i64) -> bool {
        self.setup.base.levels[l].model.reactions[j]
            .need
            .iter()
            .any(|(s, c)| self.buf[*s] / c < n_critical)
    }

    /// Mean and variance of the changes of reaction `j` (propensity `a`) in
    /// voxel `v` of level `l`, spread over the populations it changes.
    fn add_reaction_moments(&self, tau: &mut MultigridTau, l: usize, v: usize, j: usize, a: f64) {
        let setup = self.setup;
        let lay = self.lay();
        let r = &setup.base.levels[l].model.reactions[j];
        if self.same_level(l, j) {
            for &(s, d) in &r.delta {
                let (d, i) = (d as f64, tau.pop_off[s] + v);
                tau.mu[i] += d * a;
                tau.sigma2[i] += d * d * a;
            }
            return;
        }
        let add = |tau: &mut MultigridTau, s: usize, w: usize, d: f64, p: f64| {
            let i = tau.pop_off[s] + w;
            tau.mu[i] += d * a * p;
            tau.sigma2[i] += d * d * a * p;
        };
        // Where the products land: the voxels of the finest consumed species.
        tau.spots.clear();
        let finest = r
            .delta
            .iter()
            .filter(|(_, d)| *d < 0)
            .map(|(s, _)| *s)
            .min_by_key(|s| lay.level[*s])
            .filter(|s| lay.level[*s] < l);
        match finest {
            Some(f) => {
                let lf = lay.level[f];
                let total = self.agg[f][l - lf][v] as f64;
                for w in self.inside(v, l, lf) {
                    let n = self.agg[f][0][w];
                    if n > 0 {
                        tau.spots.push((lf, w, n as f64 / total));
                    }
                }
            }
            None => tau.spots.push((l, v, 1.0)),
        }
        for &(s, d) in &r.delta {
            let (ls, d) = (lay.level[s], d as f64);
            if d < 0.0 {
                let total = self.agg[s][l - ls][v] as f64;
                for w in self.inside(v, l, ls) {
                    let n = self.agg[s][0][w];
                    if n > 0 {
                        add(tau, s, w, d, n as f64 / total);
                    }
                }
                continue;
            }
            for k in 0..tau.spots.len() {
                let (lf, w, p) = tau.spots[k];
                if ls >= lf {
                    add(tau, s, setup.ancestor(w, lf, ls), d, p);
                } else {
                    let m = 8f64.powi((lf - ls) as i32);
                    for u in self.inside(w, lf, ls) {
                        add(tau, s, u, d, p / m);
                    }
                }
            }
        }
    }

    /// Rates, critical channels and the leap bound for the current counts.
    fn plan_leap(&mut self, tau: &mut MultigridTau) -> Plan {
        let setup = self.setup;
        let lay = self.lay();
        let n_sp = setup.base.n_sp;
        tau.pop_off.clear();
        let mut n_pop = 0;
        for counts in &self.agg {
            tau.pop_off.push(n_pop);
            n_pop += counts[0].len();
        }
        let n_sites = *setup.offset.last().unwrap();
        for x in [&mut tau.mu, &mut tau.sigma2, &mut tau.critical] {
            x.clear();
        }
        tau.mu.resize(n_pop, 0.0);
        tau.sigma2.resize(n_pop, 0.0);
        tau.critical.resize(n_sites, 0.0);
        let (mut a0, mut a_critical) = (0.0, 0.0);
        for i in 0..n_sites {
            let (l, v) = setup.site_of(i);
            let model = &setup.base.levels[l].model;
            let mut critical = 0.0;
            if !lay.by_level[l].is_empty() {
                self.load(l, v);
                for &j in &lay.by_level[l] {
                    let a = model.propensity(j, &self.buf);
                    if a <= 0.0 {
                        continue;
                    }
                    a0 += a;
                    if self.is_critical(l, j, tau.n_critical) {
                        critical += a;
                    } else {
                        self.add_reaction_moments(tau, l, v, j, a);
                    }
                }
            }
            let lv = &setup.base.levels[l];
            for &s in &lay.hoppers[l] {
                let count = self.agg[s][0][v];
                if count <= 0 {
                    continue;
                }
                let d = setup.base.diffusion[s];
                let is_critical = count < tau.n_critical;
                for_each_face(&lv.grid, v, |axis, w| {
                    let a = count as f64 * d * lv.inv_h2[axis];
                    a0 += a;
                    if is_critical {
                        critical += a;
                    } else {
                        let (iv, iw) = (tau.pop_off[s] + v, tau.pop_off[s] + w);
                        tau.mu[iv] -= a;
                        tau.sigma2[iv] += a;
                        tau.mu[iw] += a;
                        tau.sigma2[iw] += a;
                    }
                });
            }
            tau.critical[i] = critical;
            a_critical += critical;
        }
        let mut bound = f64::INFINITY;
        for s in 0..n_sp {
            let order = tau.orders[s];
            if order == Order::None && setup.base.diffusion[s] <= 0.0 {
                continue;
            }
            for (w, &x) in self.agg[s][0].iter().enumerate() {
                let i = tau.pop_off[s] + w;
                let (mu, sigma2) = (tau.mu[i], tau.sigma2[i]);
                if mu == 0.0 && sigma2 <= 0.0 {
                    continue;
                }
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
    fn pick_critical(&mut self, tau: &MultigridTau, target: f64) -> Option<Fire> {
        let setup = self.setup;
        let lay = self.lay();
        let mut acc = 0.0;
        let mut site = None;
        for (i, c) in tau.critical.iter().enumerate() {
            if *c > 0.0 {
                site = Some((i, acc));
                acc += c;
                if target < acc {
                    break;
                }
            }
        }
        let (i, before) = site?;
        let (l, v) = setup.site_of(i);
        let model = &setup.base.levels[l].model;
        let mut acc = before;
        let mut last = None;
        if !lay.by_level[l].is_empty() {
            self.load(l, v);
            for &j in &lay.by_level[l] {
                let a = model.propensity(j, &self.buf);
                if a > 0.0 && self.is_critical(l, j, tau.n_critical) {
                    acc += a;
                    last = Some(Fire::Reaction { l, v, j });
                    if target < acc {
                        return last;
                    }
                }
            }
        }
        let lv = &setup.base.levels[l];
        for &s in &lay.hoppers[l] {
            let count = self.agg[s][0][v];
            if count <= 0 || count >= tau.n_critical {
                continue;
            }
            let d = setup.base.diffusion[s];
            let mut found = None;
            for_each_face(&lv.grid, v, |axis, w| {
                if found.is_some() {
                    return;
                }
                acc += count as f64 * d * lv.inv_h2[axis];
                last = Some(Fire::Hop { s, v, w });
                if target < acc {
                    found = Some(Fire::Hop { s, v, w });
                }
            });
            if found.is_some() {
                return found;
            }
        }
        last
    }

    /// Apply `k` firings of reaction `j` in voxel `v` of level `l` to `trial`.
    fn apply_firings(&mut self, tau: &mut MultigridTau, l: usize, v: usize, j: usize, k: i64) {
        let setup = self.setup;
        let lay = self.lay();
        let delta = &setup.base.levels[l].model.reactions[j].delta;
        if self.same_level(l, j) {
            for &(s, d) in delta {
                tau.trial[tau.pop_off[s] + v] += d * k;
            }
            return;
        }
        for _ in 0..k {
            let mut spot = (l, v);
            for &(s, d) in delta.iter().filter(|(_, d)| *d < 0) {
                for _ in 0..-d {
                    let w = self.pick_inside(s, l, v);
                    tau.trial[tau.pop_off[s] + w] -= 1;
                    if lay.level[s] < spot.0 {
                        spot = (lay.level[s], w);
                    }
                }
            }
            for &(s, d) in delta.iter().filter(|(_, d)| *d > 0) {
                let ls = lay.level[s];
                for _ in 0..d {
                    let w = if ls >= spot.0 {
                        setup.ancestor(spot.1, spot.0, ls)
                    } else {
                        self.uniform_inside(spot.1, spot.0, ls)
                    };
                    tau.trial[tau.pop_off[s] + w] += 1;
                }
            }
        }
    }

    /// Draw the firings of a step `dt` into `tau.trial`. Returns the
    /// reactions and hops fired, or `None` if a count would go negative.
    fn draw_leap(
        &mut self,
        tau: &mut MultigridTau,
        dt: f64,
        fire: Option<Fire>,
    ) -> Option<(u64, u64)> {
        let setup = self.setup;
        let lay = self.lay();
        tau.trial.clear();
        for counts in &self.agg {
            tau.trial.extend_from_slice(&counts[0]);
        }
        let (mut reactions, mut hops) = (0u64, 0u64);
        let n_sites = *setup.offset.last().unwrap();
        for i in 0..n_sites {
            let (l, v) = setup.site_of(i);
            let model = &setup.base.levels[l].model;
            if !lay.by_level[l].is_empty() {
                self.load(l, v);
            }
            for &j in &lay.by_level[l] {
                let a = model.propensity(j, &self.buf);
                if a <= 0.0 || self.is_critical(l, j, tau.n_critical) {
                    continue;
                }
                let k = self.rng.poisson(a * dt);
                if k > 0 {
                    reactions += k as u64;
                    self.apply_firings(tau, l, v, j, k);
                }
            }
            let lv = &setup.base.levels[l];
            for &s in &lay.hoppers[l] {
                let count = self.agg[s][0][v];
                if count < tau.n_critical {
                    continue;
                }
                let d = setup.base.diffusion[s];
                let off = tau.pop_off[s];
                let rng = &mut self.rng;
                let trial = &mut tau.trial;
                for_each_face(&lv.grid, v, |axis, w| {
                    let k = rng.poisson(count as f64 * d * lv.inv_h2[axis] * dt);
                    if k > 0 {
                        hops += k as u64;
                        trial[off + v] -= k;
                        trial[off + w] += k;
                    }
                });
            }
        }
        match fire {
            Some(Fire::Reaction { l, v, j }) => {
                reactions += 1;
                self.apply_firings(tau, l, v, j, 1);
            }
            Some(Fire::Hop { s, v, w }) => {
                hops += 1;
                tau.trial[tau.pop_off[s] + v] -= 1;
                tau.trial[tau.pop_off[s] + w] += 1;
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
        let n_sites = *self.setup.offset.last().unwrap() as u64;
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
                    self.t = t_stop;
                }
                continue;
            }
            let plan = self.plan_leap(&mut tau);
            if plan.a0 <= 0.0 {
                self.t = t_stop;
                break;
            }
            let remaining = t_stop - self.t;
            let mut dt = plan.tau.min(remaining);
            if plan.a_critical > 0.0 {
                dt = dt.min(1.0 / plan.a_critical);
            }
            let expected = (plan.a0 - plan.a_critical) * dt;
            if expected < (10 + 4 * n_sites) as f64 {
                tau.exact_batch = (tau.exact_batch * 2).clamp(8 * n_sites, 512 * n_sites);
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
                let mut k = 0;
                for counts in self.agg.iter_mut() {
                    let n = counts[0].len();
                    counts[0].copy_from_slice(&tau.trial[k..k + n]);
                    k += n;
                }
                self.rebuild_aggregates();
                self.queue_stale = true;
                self.t = if to_stop { t_stop } else { self.t + dt };
                self.n_reactions += reactions;
                self.n_hops += hops;
                tau.n_leaps += 1;
                self.count_events(reactions + hops);
                break;
            }
        }
        self.tau = Some(tau);
    }
}
