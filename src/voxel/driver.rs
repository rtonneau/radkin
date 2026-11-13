//! Strang-split reaction-diffusion driver with step-doubling error control
//! (ADR 0002).
//!
//! One splitting step of `dt` is `R(dt/2) D(dt) R(dt/2)`. Each attempt
//! compares one step of `dt` with two steps of `dt/2`; the scaled max-norm
//! difference is the error estimate, the two-half-step state is kept, and dt
//! is adapted with the local order 3 of Strang splitting.
//! [`VoxelRun::advance_to`] lands exactly on the requested time, so injections
//! and snapshots happen between calls.

use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::solver::SolveOptions;
use crate::voxel::diffusion::{Diffusion, Scheme};
use crate::voxel::grid::VoxelState;
use crate::voxel::reaction::reaction_step;

/// Step-size control of the splitting.
#[derive(Debug, Clone, PartialEq)]
pub struct StepControl {
    pub rtol: f64,
    /// Absolute tolerance [µmol/L].
    pub atol: f64,
    /// First splitting step [s].
    pub dt_initial: f64,
    pub dt_min: f64,
    pub dt_max: f64,
    /// Largest growth / smallest shrink factor per step.
    pub grow: f64,
    pub shrink: f64,
    pub safety: f64,
    /// Backward Euler substeps of the diffusion after an injection (Rannacher).
    pub rannacher_steps: usize,
}

impl Default for StepControl {
    fn default() -> Self {
        Self {
            rtol: 1e-3,
            atol: 1e-6,
            dt_initial: 1e-12,
            dt_min: 1e-18,
            dt_max: f64::INFINITY,
            grow: 2.0,
            shrink: 0.2,
            safety: 0.9,
            rannacher_steps: 4,
        }
    }
}

/// One accepted splitting step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepRecord {
    pub t: f64,
    pub dt: f64,
    pub error: f64,
}

/// A voxel reaction-diffusion run.
#[derive(Debug, Clone)]
pub struct VoxelRun {
    pub mat: SimMatrices,
    pub diffusion: Diffusion,
    pub state: VoxelState,
    pub t: f64,
    pub control: StepControl,
    /// Per-voxel rsoda options (`t_eval` ignored).
    pub solve: SolveOptions,
    pub steps: Vec<StepRecord>,
    pub rejected: usize,
    /// When `Some`, the domain sum of every species [µmol/L summed over
    /// voxels] after each accepted step, as `(t, sums)`.
    pub totals: Option<Vec<(f64, Vec<f64>)>>,
    dt: f64,
    rannacher_pending: bool,
}

impl VoxelRun {
    pub fn new(
        mat: SimMatrices,
        diffusion: Diffusion,
        state: VoxelState,
        t0: f64,
        control: StepControl,
        solve: SolveOptions,
    ) -> Self {
        let dt = control.dt_initial;
        Self {
            mat,
            diffusion,
            state,
            t: t0,
            control,
            solve,
            steps: Vec::new(),
            rejected: 0,
            totals: None,
            dt,
            rannacher_pending: true,
        }
    }

    /// Add `amount` [µmol/L, same layout as the state] and damp the next
    /// diffusion step (Rannacher). dt restarts from `dt_initial`.
    pub fn inject(&mut self, amount: &VoxelState) {
        assert_eq!(amount.data.len(), self.state.data.len());
        for (s, a) in self.state.data.iter_mut().zip(&amount.data) {
            *s += a;
        }
        self.rannacher_pending = true;
        self.dt = self.control.dt_initial;
        self.record_totals();
    }

    /// Start recording domain sums (records the current state first).
    pub fn track_totals(&mut self) {
        self.totals = Some(Vec::new());
        self.record_totals();
    }

    /// Domain sum of every species [µmol/L summed over voxels].
    pub fn domain_sums(&self) -> Vec<f64> {
        (0..self.state.n_species)
            .map(|s| self.state.field(s).iter().sum())
            .collect()
    }

    fn record_totals(&mut self) {
        if self.totals.is_some() {
            let sums = self.domain_sums();
            if let Some(v) = self.totals.as_mut() {
                v.push((self.t, sums));
            }
        }
    }

    fn strang(&mut self, state: &mut VoxelState, t: f64, dt: f64, damped: bool) -> Result<()> {
        reaction_step(&self.mat, state, t, dt / 2.0, &self.solve)?;
        if damped {
            self.diffusion
                .step_rannacher(state, dt, self.control.rannacher_steps);
        } else {
            self.diffusion.step(state, dt, Scheme::CrankNicolson);
        }
        reaction_step(&self.mat, state, t + dt / 2.0, dt / 2.0, &self.solve)
    }

    fn error_norm(&self, big: &VoxelState, small: &VoxelState) -> f64 {
        let c = &self.control;
        big.data
            .iter()
            .zip(&small.data)
            .map(|(a, b)| (a - b).abs() / (c.atol + c.rtol * b.abs()))
            .fold(0.0, f64::max)
    }

    /// Advance to exactly `t_target`.
    pub fn advance_to(&mut self, t_target: f64) -> Result<()> {
        while self.t < t_target {
            let remaining = t_target - self.t;
            let mut dt = self.dt.min(self.control.dt_max);
            let landing = dt >= remaining * (1.0 - 1e-12);
            if landing {
                dt = remaining;
            }
            let damped = self.rannacher_pending;
            let mut big = self.state.clone();
            self.strang(&mut big, self.t, dt, damped)?;
            let mut small = self.state.clone();
            self.strang(&mut small, self.t, dt / 2.0, damped)?;
            self.strang(&mut small, self.t + dt / 2.0, dt / 2.0, false)?;
            let err = self.error_norm(&big, &small);
            let c = &self.control;
            let factor = if err == 0.0 {
                c.grow
            } else {
                (c.safety * err.powf(-1.0 / 3.0)).clamp(c.shrink, c.grow)
            };
            if err <= 1.0 {
                self.state = small;
                self.t = if landing { t_target } else { self.t + dt };
                self.steps.push(StepRecord {
                    t: self.t,
                    dt,
                    error: err,
                });
                self.rannacher_pending = false;
                self.record_totals();
                // A landing step was cut short: keep the earlier proposal if larger.
                let next = dt * factor;
                self.dt = if landing { self.dt.max(next) } else { next };
            } else {
                self.rejected += 1;
                self.dt = dt * factor;
                if self.dt < self.control.dt_min {
                    return Err(RadkinError::Value(format!(
                        "voxel splitting step fell below dt_min={:e} s at t={:e} s (error {err:e})",
                        self.control.dt_min, self.t
                    )));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::solver::{DoseRateSource, solve};
    use crate::voxel::grid::VoxelGrid;
    use crate::voxel::reaction::tests::{opts, water_env};

    #[test]
    fn uniform_seed_matches_unit_cell() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let t_end = 1e-3;
        let o = SolveOptions {
            t_eval: Some(vec![t_end]),
            ..opts()
        };
        let reference = solve(&mat, &DoseRateSource::None, (0.0, t_end), &y0, &o)
            .unwrap()
            .y[0]
            .clone();
        let grid = VoxelGrid::new([3, 2, 2], [5e-9; 3], [0.0; 3]).unwrap();
        let coeffs = vec![2e-9; y0.len()];
        let st = VoxelState::uniform(&y0, grid.n_voxels());
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid, coeffs),
            st,
            0.0,
            StepControl::default(),
            opts(),
        );
        run.advance_to(t_end / 2.0).unwrap();
        run.advance_to(t_end).unwrap();
        assert_eq!(run.t, t_end);
        let mut y = vec![0.0; y0.len()];
        for v in 0..run.state.n_voxels {
            run.state.gather(v, &mut y);
            for (a, b) in y.iter().zip(&reference) {
                let tol = 1e-6 * b.abs().max(1e-6);
                assert!((a - b).abs() <= tol, "voxel {v}: {a} vs {b}");
            }
        }
        assert!(run.steps.len() < 200, "{} steps", run.steps.len());
    }

    fn ab_run(rtol: f64) -> Vec<f64> {
        let cfg = Config::from_toml_str(
            "[parameters]\n[concentrations.initial]\nA = 0.0\nB = 0.0\n[[reactions.k_reaction]]\nreaction = \"A + B -> C\"\nk_value = 1e10\n",
            None,
        )
        .unwrap();
        let env = SimEnv::from_config(&cfg).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let ia = env.species.index_of("A").unwrap();
        let ib = env.species.index_of("B").unwrap();
        let n = env.species.n_dyn();
        let grid = VoxelGrid::new([20, 1, 1], [5e-9; 3], [0.0; 3]).unwrap();
        let mut st = VoxelState::zeros(n, grid.n_voxels());
        for v in 0..10 {
            st.field_mut(ia)[v] = 1e3; // 1 mM
            st.field_mut(ib)[v + 10] = 1e3;
        }
        let mut coeffs = vec![0.0; n];
        coeffs[ia] = 2e-9;
        coeffs[ib] = 1e-9;
        let control = StepControl {
            rtol,
            atol: 1e-3,
            ..StepControl::default()
        };
        let mut run = VoxelRun::new(mat, Diffusion::new(grid, coeffs), st, 0.0, control, opts());
        run.advance_to(2e-7).unwrap();
        run.state.data
    }

    #[test]
    fn error_decreases_with_tolerance() {
        let reference = ab_run(1e-7);
        let err = |rtol| {
            ab_run(rtol)
                .iter()
                .zip(&reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        };
        let (coarse, fine) = (err(1e-2), err(1e-4));
        assert!(fine < coarse, "fine {fine} coarse {coarse}");
        assert!(fine < 1.0, "fine {fine} µM");
    }

    #[test]
    fn injection_is_added_and_conserved_by_diffusion() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let n = env.species.n_dyn();
        let grid = VoxelGrid::new([5, 5, 1], [5e-9; 3], [0.0; 3]).unwrap();
        let st = VoxelState::zeros(n, grid.n_voxels());
        let mut coeffs = vec![0.0; n];
        let ih2 = env.species.index_of("H2").unwrap();
        coeffs[ih2] = 5e-9;
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid.clone(), coeffs),
            st,
            0.0,
            StepControl::default(),
            opts(),
        );
        run.advance_to(1e-12).unwrap();
        let mut add = VoxelState::zeros(n, grid.n_voxels());
        add.field_mut(ih2)[grid.index(2, 2, 0)] = 100.0;
        run.inject(&add);
        run.advance_to(1e-9).unwrap();
        let total: f64 = run.state.field(ih2).iter().sum();
        // H2 is nearly inert at these concentrations over 1 ns.
        assert!((total - 100.0).abs() < 1e-3, "{total}");
        assert!(run.state.field(ih2).iter().all(|v| *v >= -1e-9));
        assert!(run.state.field(ih2)[grid.index(0, 0, 0)] > 0.0);
    }
}
