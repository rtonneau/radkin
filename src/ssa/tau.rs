//! Tau-leaping of one well-mixed volume (opt-in, `ssa.tau_leap = true`).
//!
//! Explicit Poisson tau-leaping with the step selection of Cao, Gillespie
//! and Petzold (J. Chem. Phys. 124, 044109, 2006) and their handling of
//! critical reactions (J. Chem. Phys. 123, 054104, 2005):
//!
//! 1. A reaction is critical when it has a positive propensity and fewer
//!    than [`TauLeap::n_critical`] firings left before it exhausts a species
//!    it consumes. Critical reactions fire one at a time, as in the exact SSA.
//! 2. The leap `tau'` bounds the relative change of every propensity by
//!    [`TauLeap::epsilon`], from the mean and variance of the change of each
//!    reactant species over the non-critical reactions.
//! 3. If `tau'` is below [`TauLeap::ssa_threshold`] expected events, leaping
//!    does not pay: [`TauLeap::ssa_steps`] exact SSA steps are taken instead.
//! 4. Otherwise the time `tau''` of the next critical reaction is drawn, the
//!    step is the smaller of the two, every non-critical reaction fires a
//!    Poisson number of times, and one critical reaction fires if `tau''`
//!    was the smaller.
//! 5. A step that would make a count negative is rejected and redone with
//!    `tau'` halved.
//!
//! A leap never crosses the end of the span it is asked to advance over, so
//! beam transitions and output times are honoured exactly.

use crate::ssa::exact::{Cell, Stepper};
use crate::ssa::model::{Law, Slot, SsaModel};
use crate::ssa::rng::SsaRng;

/// How a species enters the highest-order reactions it is a reactant of
/// (the `g_i` of Cao, Gillespie and Petzold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Order {
    /// Not a reactant: its count bounds no propensity.
    None,
    First,
    Second,
    /// Second order with two molecules of the species.
    SecondSelf,
    Third,
    /// Third order with two molecules of the species.
    ThirdTwo,
    /// Third order with three molecules of the species.
    ThirdSelf,
    /// Higher order: `g = 2 * order` (conservative).
    High(u32),
}

impl Order {
    fn rank(&self) -> (u32, u32) {
        match self {
            Order::None => (0, 0),
            Order::First => (1, 1),
            Order::Second => (2, 1),
            Order::SecondSelf => (2, 2),
            Order::Third => (3, 1),
            Order::ThirdTwo => (3, 2),
            Order::ThirdSelf => (3, 3),
            Order::High(n) => (*n, 4),
        }
    }

    fn of(order: u32, multiplicity: u32) -> Self {
        match (order, multiplicity) {
            (0 | 1, _) => Order::First,
            (2, 1) => Order::Second,
            (2, _) => Order::SecondSelf,
            (3, 1) => Order::Third,
            (3, 2) => Order::ThirdTwo,
            (3, _) => Order::ThirdSelf,
            (n, _) => Order::High(n),
        }
    }

    /// `g_i` for a count `x`.
    pub(crate) fn g(&self, x: i64) -> f64 {
        let x = x as f64;
        match self {
            Order::None | Order::First => 1.0,
            Order::Second => 2.0,
            Order::SecondSelf if x >= 2.0 => 2.0 + 1.0 / (x - 1.0),
            Order::SecondSelf => 3.0,
            Order::Third => 3.0,
            Order::ThirdTwo if x >= 2.0 => 1.5 * (2.0 + 1.0 / (x - 1.0)),
            Order::ThirdTwo => 4.5,
            Order::ThirdSelf if x >= 3.0 => 3.0 + 1.0 / (x - 1.0) + 2.0 / (x - 2.0),
            Order::ThirdSelf => 5.5,
            Order::High(n) => 2.0 * *n as f64,
        }
    }
}

/// Highest-order role of every species of `model`.
pub(crate) fn species_orders(model: &SsaModel) -> Vec<Order> {
    let mut orders = vec![Order::None; model.n_dyn];
    let mut raise = |s: usize, o: Order| {
        if o.rank() > orders[s].rank() {
            orders[s] = o;
        }
    };
    for r in &model.reactions {
        let order = r.order.max(1);
        match &r.law {
            Law::MassAction { reactants, .. } => {
                for &(s, mult) in reactants {
                    raise(s, Order::of(order, mult));
                }
            }
            Law::Michaelis {
                enzyme, substrate, ..
            } => {
                for slot in [enzyme, substrate] {
                    if let Slot::Dyn(s) = slot {
                        raise(*s, Order::of(order, 1));
                    }
                }
            }
        }
        // A species only consumed (a "negative" product) switches the
        // reaction off when it runs out: it bounds the step as a reactant.
        for &(s, _) in &r.need {
            raise(s, Order::of(order, 1));
        }
    }
    orders
}

/// Tau-leaping stepper of one trajectory.
#[derive(Debug, Clone)]
pub struct TauLeap {
    /// Bound on the relative change of a propensity over a leap.
    pub epsilon: f64,
    /// A reaction with fewer firings left than this is critical.
    pub n_critical: i64,
    /// A leap expected to fire fewer events than this is replaced by exact steps.
    pub ssa_threshold: f64,
    /// Exact steps taken when leaping does not pay.
    pub ssa_steps: usize,
    /// Leaps taken.
    pub n_leaps: u64,
    /// Leaps rejected for a negative count.
    pub n_rejected: u64,
    /// Exact steps taken.
    pub n_exact: u64,
    orders: Vec<Order>,
    critical: Vec<bool>,
    fired: Vec<i64>,
    fired_sources: Vec<i64>,
    mu: Vec<f64>,
    sigma2: Vec<f64>,
    trial: Vec<i64>,
}

impl TauLeap {
    /// Stepper with the default epsilon (0.03).
    pub fn new(model: &SsaModel) -> Self {
        Self::with_epsilon(model, crate::ssa::config::DEFAULT_TAU_EPSILON)
    }

    /// Stepper bounding the relative change of a propensity by `epsilon`.
    pub fn with_epsilon(model: &SsaModel, epsilon: f64) -> Self {
        Self {
            epsilon,
            n_critical: 10,
            ssa_threshold: 10.0,
            ssa_steps: 100,
            n_leaps: 0,
            n_rejected: 0,
            n_exact: 0,
            orders: species_orders(model),
            critical: vec![false; model.reactions.len()],
            fired: vec![0; model.reactions.len()],
            fired_sources: vec![0; model.sources.len()],
            mu: vec![0.0; model.n_dyn],
            sigma2: vec![0.0; model.n_dyn],
            trial: vec![0; model.n_dyn],
        }
    }

    /// Mark the critical reactions; returns their total propensity.
    fn mark_critical(&mut self, cell: &Cell<'_>) -> f64 {
        let mut total = 0.0;
        for (j, r) in cell.model.reactions.iter().enumerate() {
            let a = cell.a[j];
            let mut critical = false;
            if a > 0.0 {
                for &(s, c) in &r.need {
                    if cell.counts[s] / c < self.n_critical {
                        critical = true;
                        break;
                    }
                }
            }
            self.critical[j] = critical;
            if critical {
                total += a;
            }
        }
        total
    }

    /// The leap bounding the change of every propensity (`tau'`).
    fn leap_bound(&mut self, cell: &Cell<'_>, dose_rate: f64) -> f64 {
        self.mu.fill(0.0);
        self.sigma2.fill(0.0);
        for (j, r) in cell.model.reactions.iter().enumerate() {
            let a = cell.a[j];
            if a <= 0.0 || self.critical[j] {
                continue;
            }
            for &(s, d) in &r.delta {
                let d = d as f64;
                self.mu[s] += d * a;
                self.sigma2[s] += d * d * a;
            }
        }
        if dose_rate > 0.0 {
            for src in &cell.model.sources {
                let a = src.per_gray * dose_rate;
                self.mu[src.species] += a;
                self.sigma2[src.species] += a;
            }
        }
        let mut tau = f64::INFINITY;
        for (s, order) in self.orders.iter().enumerate() {
            if *order == Order::None {
                continue;
            }
            let x = cell.counts[s];
            let bound = (self.epsilon * x as f64 / order.g(x)).max(1.0);
            if self.mu[s] != 0.0 {
                tau = tau.min(bound / self.mu[s].abs());
            }
            if self.sigma2[s] > 0.0 {
                tau = tau.min(bound * bound / self.sigma2[s]);
            }
        }
        tau
    }

    /// Draw the firings of a step `tau` into `trial`; false if a count
    /// would go negative.
    fn draw(
        &mut self,
        cell: &Cell<'_>,
        tau: f64,
        dose_rate: f64,
        critical_fire: Option<usize>,
        rng: &mut SsaRng,
    ) -> bool {
        let model = cell.model;
        self.trial.copy_from_slice(&cell.counts);
        for (j, r) in model.reactions.iter().enumerate() {
            let a = cell.a[j];
            let mut k = if a > 0.0 && !self.critical[j] {
                rng.poisson(a * tau)
            } else {
                0
            };
            if critical_fire == Some(j) {
                k += 1;
            }
            self.fired[j] = k;
            if k > 0 {
                for &(s, d) in &r.delta {
                    self.trial[s] += d * k;
                }
            }
        }
        for (i, src) in model.sources.iter().enumerate() {
            let k = if dose_rate > 0.0 {
                rng.poisson(src.per_gray * dose_rate * tau)
            } else {
                0
            };
            self.fired_sources[i] = k;
            self.trial[src.species] += k;
        }
        self.trial.iter().all(|n| *n >= 0)
    }
}

impl Stepper for TauLeap {
    fn advance(&mut self, cell: &mut Cell<'_>, t_stop: f64, dose_rate: f64, rng: &mut SsaRng) {
        while cell.t < t_stop {
            let a_reactions: f64 = cell.a.iter().sum();
            let a_sources: f64 = if dose_rate > 0.0 {
                let per_gray: f64 = cell.model.sources.iter().map(|s| s.per_gray).sum();
                per_gray * dose_rate
            } else {
                0.0
            };
            let a0 = a_reactions + a_sources;
            if a0 <= 0.0 {
                cell.t = t_stop;
                return;
            }
            let a_critical = self.mark_critical(cell);
            let mut tau_leap = self.leap_bound(cell, dose_rate);
            if tau_leap < self.ssa_threshold / a0 {
                for _ in 0..self.ssa_steps {
                    if !cell.step(t_stop, dose_rate, rng) {
                        break;
                    }
                    self.n_exact += 1;
                }
                continue;
            }
            loop {
                let remaining = t_stop - cell.t;
                let tau_critical = if a_critical > 0.0 {
                    rng.exp1() / a_critical
                } else {
                    f64::INFINITY
                };
                // The smallest of the leap, the next critical firing and the
                // end of the span decides the step.
                let (tau, fire, to_stop) = if remaining <= tau_leap && remaining <= tau_critical {
                    (remaining, false, true)
                } else if tau_critical < tau_leap {
                    (tau_critical, true, false)
                } else {
                    (tau_leap, false, false)
                };
                let critical_fire = fire.then(|| {
                    let target = rng.uniform() * a_critical;
                    let mut acc = 0.0;
                    let mut last = 0;
                    for j in 0..cell.a.len() {
                        if self.critical[j] {
                            acc += cell.a[j];
                            last = j;
                            if target < acc {
                                break;
                            }
                        }
                    }
                    last
                });
                if !self.draw(cell, tau, dose_rate, critical_fire, rng) {
                    self.n_rejected += 1;
                    tau_leap *= 0.5;
                    continue;
                }
                cell.counts.copy_from_slice(&self.trial);
                cell.t = if to_stop { t_stop } else { cell.t + tau };
                let fired: i64 =
                    self.fired.iter().sum::<i64>() + self.fired_sources.iter().sum::<i64>();
                cell.n_events += fired as u64;
                cell.refresh();
                self.n_leaps += 1;
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::matrices::SimMatrices;

    fn model(toml: &str, volume: f64) -> (SimEnv, SsaModel) {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        (env, SsaModel::new(&mat, volume).unwrap())
    }

    const CHEMISTRY: &str = r#"
        [[reactions.k_reaction]]
        reaction = "H_r -> H2"
        k_value = 3.0e3

        [[reactions.k_reaction]]
        reaction = "2 OH_r -> H2O2"
        k_value = 5.5e9

        [[reactions.k_reaction]]
        reaction = "e_aq + OH_r -> OH_minus"
        k_value = 3.0e10
    "#;

    #[test]
    fn orders_follow_the_highest_order_reaction() {
        let (env, m) = model(CHEMISTRY, 1e-15);
        let orders = species_orders(&m);
        let of = |l: &str| orders[env.species.index_of(l).unwrap()];
        assert_eq!(of("H_r"), Order::First);
        assert_eq!(of("OH_r"), Order::SecondSelf);
        assert_eq!(of("e_aq"), Order::Second);
        assert_eq!(of("H2"), Order::None);
        assert_eq!(of("H2O2"), Order::None);
        assert_eq!(Order::SecondSelf.g(101), 2.01);
        assert_eq!(Order::First.g(5), 1.0);
        assert!((Order::ThirdSelf.g(4) - (3.0 + 1.0 / 3.0 + 1.0)).abs() < 1e-15);
    }

    #[test]
    fn counts_stay_non_negative_and_conserved() {
        let (env, m) = model(CHEMISTRY, 1e-15);
        let idx = |l: &str| env.species.index_of(l).unwrap();
        let (h, h2, oh, h2o2, e, ohm) = (
            idx("H_r"),
            idx("H2"),
            idx("OH_r"),
            idx("H2O2"),
            idx("e_aq"),
            idx("OH_minus"),
        );
        for stream in 0..50 {
            let mut n0 = vec![0i64; m.n_dyn];
            n0[h] = 30_000;
            n0[oh] = 50_000;
            n0[e] = 20_000;
            let mut cell = Cell::new(&m, n0, 0.0);
            let mut tau = TauLeap::new(&m);
            let mut rng = SsaRng::new(9, stream);
            let mut t = 1e-9;
            while t < 1e-2 {
                tau.advance(&mut cell, t, 0.0, &mut rng);
                assert_eq!(cell.t, t, "the step ends exactly at the requested time");
                assert!(cell.counts.iter().all(|n| *n >= 0));
                assert_eq!(cell.counts[h] + cell.counts[h2], 30_000);
                assert_eq!(cell.counts[e] + cell.counts[ohm], 20_000);
                assert_eq!(
                    cell.counts[oh] + 2 * cell.counts[h2o2] + cell.counts[ohm],
                    50_000
                );
                t *= 3.0;
            }
            assert!(tau.n_leaps > 0, "leaps were taken");
            // Far fewer steps than events.
            assert!(tau.n_leaps + tau.n_exact < cell.n_events / 10);
        }
    }

    #[test]
    fn small_counts_fall_back_to_exact_steps() {
        let (env, m) = model(CHEMISTRY, 1e-15);
        let mut n0 = vec![0i64; m.n_dyn];
        n0[env.species.index_of("H_r").unwrap()] = 8;
        let mut cell = Cell::new(&m, n0, 0.0);
        let mut tau = TauLeap::new(&m);
        let mut rng = SsaRng::new(1, 0);
        tau.advance(&mut cell, 1.0, 0.0, &mut rng);
        assert_eq!(cell.counts[env.species.index_of("H2").unwrap()], 8);
        assert_eq!(cell.n_events, 8, "every decay is one event");
        assert_eq!(tau.n_rejected, 0);
    }
}
