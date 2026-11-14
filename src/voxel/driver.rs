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
use crate::voxel::reaction::{integrate_one, reaction_step_background};

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
    /// With `totals` and a background: the background vector times the
    /// number of voxels at the same records (what the domain would hold
    /// without seeds).
    pub background_totals: Option<Vec<Vec<f64>>>,
    /// Background vector [µmol/L]: the uniform solution without seeds,
    /// advanced with the state. With `skip_background`, voxels that match it
    /// share its integration (see [`reaction_step_background`]).
    pub background: Option<Vec<f64>>,
    /// Advance background voxels with the shared integration (default true).
    pub skip_background: bool,
    /// Background voxels summed over all reaction half-steps (rejected
    /// attempts included).
    pub background_voxel_steps: u64,
    /// Voxels summed over all reaction half-steps.
    pub voxel_steps: u64,
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
            background_totals: None,
            background: None,
            skip_background: true,
            background_voxel_steps: 0,
            voxel_steps: 0,
            dt,
            rannacher_pending: true,
        }
    }

    /// Splitting step size proposed for the next attempt [s].
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Set the next splitting step size [s] (e.g. when resuming a checkpoint).
    pub fn set_dt(&mut self, dt: f64) {
        self.dt = dt;
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
        self.background_totals = self.background.as_ref().map(|_| Vec::new());
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
            if let (Some(v), Some(bg)) = (self.background_totals.as_mut(), &self.background) {
                let n = self.state.n_voxels as f64;
                v.push(bg.iter().map(|b| b * n).collect());
            }
        }
    }

    fn strang(
        &mut self,
        state: &mut VoxelState,
        bg: &mut Option<Vec<f64>>,
        t: f64,
        dt: f64,
        damped: bool,
    ) -> Result<()> {
        self.half_reaction(state, bg, t, dt / 2.0)?;
        let shift = if self.skip_background {
            bg.as_deref()
        } else {
            None
        };
        self.diffuse(state, shift, dt, damped);
        self.half_reaction(state, bg, t + dt / 2.0, dt / 2.0)
    }

    /// Diffusion over `dt`. With a background, the deviation from it is
    /// diffused (zero-flux diffusion leaves a uniform field unchanged), so
    /// background voxels stay exactly background instead of drifting by
    /// round-off.
    fn diffuse(&mut self, state: &mut VoxelState, bg: Option<&[f64]>, dt: f64, damped: bool) {
        let shift = |state: &mut VoxelState, sign: f64| {
            if let Some(bg) = bg {
                for (s, b) in bg.iter().enumerate() {
                    for x in state.field_mut(s) {
                        *x += sign * b;
                    }
                }
            }
        };
        shift(state, -1.0);
        if damped {
            self.diffusion
                .step_rannacher(state, dt, self.control.rannacher_steps);
        } else {
            self.diffusion.step(state, dt, Scheme::CrankNicolson);
        }
        shift(state, 1.0);
    }

    fn half_reaction(
        &mut self,
        state: &mut VoxelState,
        bg: &mut Option<Vec<f64>>,
        t: f64,
        h: f64,
    ) -> Result<()> {
        let n_bg = if self.skip_background {
            reaction_step_background(&self.mat, state, bg.as_mut(), t, h, &self.solve)?
        } else {
            if let Some(b) = bg.as_mut() {
                integrate_one(&self.mat, b, t, h, &self.solve)?;
            }
            reaction_step_background(&self.mat, state, None, t, h, &self.solve)?
        };
        self.background_voxel_steps += n_bg as u64;
        self.voxel_steps += state.n_voxels as u64;
        Ok(())
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
            let mut big_bg = self.background.clone();
            self.strang(&mut big, &mut big_bg, self.t, dt, damped)?;
            let mut small = self.state.clone();
            let mut small_bg = self.background.clone();
            self.strang(&mut small, &mut small_bg, self.t, dt / 2.0, damped)?;
            self.strang(
                &mut small,
                &mut small_bg,
                self.t + dt / 2.0,
                dt / 2.0,
                false,
            )?;
            let err = self.error_norm(&big, &small);
            let c = &self.control;
            let factor = if err == 0.0 {
                c.grow
            } else {
                (c.safety * err.powf(-1.0 / 3.0)).clamp(c.shrink, c.grow)
            };
            if err <= 1.0 {
                self.state = small;
                self.background = small_bg;
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

    fn seeded_run(background: bool) -> VoxelRun {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let n = y0.len();
        let grid = VoxelGrid::new([10, 10, 3], [5e-9; 3], [0.0; 3]).unwrap();
        let coeffs = vec![3e-9; n];
        let control = StepControl {
            rtol: 1e-3,
            atol: 1e-6,
            ..StepControl::default()
        };
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid.clone(), coeffs),
            VoxelState::uniform(&y0, grid.n_voxels()),
            0.0,
            control,
            opts(),
        );
        if background {
            run.background = Some(y0.clone());
        }
        run.advance_to(1e-12).unwrap();
        let mut add = VoxelState::zeros(n, grid.n_voxels());
        for s in ["e_aq", "OH_r", "H3O_plus"] {
            let i = env.species.index_of(s).unwrap();
            add.field_mut(i)[grid.index(1, 1, 1)] = 50.0;
            add.field_mut(i)[grid.index(2, 1, 1)] = 20.0;
        }
        run.inject(&add);
        run.advance_to(2e-11).unwrap();
        run
    }

    #[test]
    fn background_voxels_match_full_integration() {
        let full = seeded_run(false);
        let skip = seeded_run(true);
        assert_eq!(full.background_voxel_steps, 0);
        assert!(
            skip.background_voxel_steps > 0,
            "no background voxel was skipped"
        );
        assert!(skip.background_voxel_steps < skip.voxel_steps);
        assert_eq!(full.steps.len(), skip.steps.len());
        for (a, b) in skip.state.data.iter().zip(&full.state.data) {
            let tol = 1e-8 * b.abs() + 1e-10;
            assert!((a - b).abs() <= tol, "{a} vs {b}");
        }
    }

    #[test]
    fn unseeded_grid_is_all_background() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let grid = VoxelGrid::new([6, 5, 4], [5e-9; 3], [0.0; 3]).unwrap();
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid.clone(), vec![2e-9; y0.len()]),
            VoxelState::uniform(&y0, grid.n_voxels()),
            0.0,
            StepControl::default(),
            opts(),
        );
        run.background = Some(y0);
        run.advance_to(1e-6).unwrap();
        assert!(run.voxel_steps > 0);
        assert_eq!(run.background_voxel_steps, run.voxel_steps);
    }

    #[test]
    fn diffusion_keeps_a_uniform_background_exact() {
        // Magnitudes of an aerated background: O2, H3O+, OH-, a trace species.
        let bg = vec![52.0, 0.1, 0.1, 1e-9, 0.0];
        let grid = VoxelGrid::new([7, 9, 11], [3e-9, 5e-9, 4e-9], [0.0; 3]).unwrap();
        let coeffs = vec![2.4e-9, 9.5e-9, 5.3e-9, 7e-9, 0.0];
        let env = water_env();
        let mut run = VoxelRun::new(
            SimMatrices::build(&env).unwrap(),
            Diffusion::new(grid.clone(), coeffs),
            VoxelState::zeros(bg.len(), grid.n_voxels()),
            0.0,
            StepControl::default(),
            opts(),
        );
        let mut st = VoxelState::uniform(&bg, grid.n_voxels());
        run.diffuse(&mut st, Some(&bg), 1e-12, true);
        for k in 0..40 {
            run.diffuse(&mut st, Some(&bg), 1e-12 * 1.7f64.powi(k), false);
        }
        assert_eq!(st.data, VoxelState::uniform(&bg, grid.n_voxels()).data);
    }
}
