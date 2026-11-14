//! Strang-split reaction-diffusion driver with step-doubling error control
//! (ADR 0002).
//!
//! One splitting step of `dt` is `R(dt/2) D(dt) R(dt/2)`. Each attempt
//! compares one step of `dt` with two steps of `dt/2`; the scaled max-norm
//! difference is the error estimate, the two-half-step state is kept, and dt
//! is adapted with the local order 3 of Strang splitting.
//! [`VoxelRun::advance_to`] lands exactly on the requested time, so injections
//! and snapshots happen between calls.
//!
//! Global coarsening ([`VoxelRun::set_coarsening`], off by default): the whole
//! grid sits at one coarsening level L, each voxel covering 2^L x 2^L x 2^L fine
//! voxels. After each accepted step, if L is below the maximum and the state
//! passes [`can_coarsen`] with the control tolerances, the grid drops one level
//! (volume average, see [`crate::voxel::coarsen`]); the dt proposal is kept.
//! [`VoxelRun::inject`] refines to level 0 (piecewise constant) before adding a
//! deposit.

use std::borrow::Cow;

use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::solver::SolveOptions;
use crate::voxel::coarsen::{can_coarsen, coarse_grid, prolonged, restricted};
use crate::voxel::diffusion::{Diffusion, Scheme};
use crate::voxel::grid::{VoxelGrid, VoxelState};
use crate::voxel::reaction::{CarriedSteps, integrate_one, reaction_step_background};

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
///
/// `state`, `diffusion` and the carried rsoda steps live on the grid of the
/// current coarsening level ([`Self::level`]); [`Self::fine_state`] gives the
/// state on the fine grid.
#[derive(Debug, Clone)]
pub struct VoxelRun {
    pub mat: SimMatrices,
    /// Diffusion operator on the current level's grid.
    pub diffusion: Diffusion,
    /// State on the current level's grid (see [`Self::fine_state`]).
    pub state: VoxelState,
    pub t: f64,
    pub control: StepControl,
    /// Per-voxel rsoda options (`t_eval` ignored).
    pub solve: SolveOptions,
    pub steps: Vec<StepRecord>,
    pub rejected: usize,
    /// When `Some`, the domain sum of every species [µmol/L summed over fine
    /// voxels, see [`Self::domain_sums`]] after each accepted step, as `(t, sums)`.
    pub totals: Option<Vec<(f64, Vec<f64>)>>,
    /// With `totals` and a background: the background vector times the
    /// number of fine voxels at the same records (what the domain would hold
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
    /// Voxels (of the level they ran on) summed over all reaction half-steps.
    pub voxel_steps: u64,
    /// `(t, level)` after every change of coarsening level, in order [s].
    pub level_history: Vec<(f64, u32)>,
    dt: f64,
    rannacher_pending: bool,
    /// rsoda step sizes carried between reaction half-steps, along the
    /// accepted trajectory (not stored in checkpoints).
    carry: CarriedSteps,
    /// Current coarsening level (0 = the fine grid).
    level: u32,
    /// Coarsening levels allowed (0 = off).
    max_coarsen_levels: u32,
    /// Tolerance factor of the coarsening test ([`can_coarsen`]).
    coarsen_tol_factor: f64,
    /// Grid of every level up to `max_coarsen_levels`; `grids[0]` is the fine grid.
    grids: Vec<VoxelGrid>,
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
        let carry = CarriedSteps::new(state.n_voxels);
        let grids = vec![diffusion.grid().clone()];
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
            level_history: Vec::new(),
            dt,
            rannacher_pending: true,
            carry,
            level: 0,
            max_coarsen_levels: 0,
            coarsen_tol_factor: 0.1,
            grids,
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

    /// Allow up to `max_levels` coarsening levels, tested with the factor
    /// `tol_factor` (`[voxel] max_coarsen_levels` and `coarsen_tol_factor`; 0
    /// levels turns coarsening off). Every fine grid dimension must be
    /// divisible by `2^max_levels`. A run above the new maximum is refined to it.
    pub fn set_coarsening(&mut self, max_levels: u32, tol_factor: f64) -> Result<()> {
        if !(tol_factor.is_finite() && tol_factor > 0.0) {
            return Err(RadkinError::Value(format!(
                "voxel coarsening: tolerance factor must be finite and > 0, got {tol_factor}"
            )));
        }
        let fine = &self.grids[0];
        let block = 1usize.checked_shl(max_levels).ok_or_else(|| {
            RadkinError::Value(format!("voxel coarsening: {max_levels} levels is too many"))
        })?;
        for (axis, n) in ["x", "y", "z"].iter().zip(fine.dims) {
            if n % block != 0 {
                return Err(RadkinError::Value(format!(
                    "voxel coarsening: grid dims along {axis} ({n}) must be divisible by \
                     2^{max_levels} = {block}"
                )));
            }
        }
        let mut grids = vec![fine.clone()];
        for l in 0..max_levels as usize {
            grids.push(coarse_grid(&grids[l])?);
        }
        self.refine_to(self.level.min(max_levels));
        self.grids = grids;
        self.max_coarsen_levels = max_levels;
        self.coarsen_tol_factor = tol_factor;
        Ok(())
    }

    /// Current coarsening level (0 = the fine grid).
    pub fn level(&self) -> u32 {
        self.level
    }

    /// The fine grid (level 0).
    pub fn fine_grid(&self) -> &VoxelGrid {
        &self.grids[0]
    }

    /// The state on the fine grid: [`Self::state`] itself at level 0, else the
    /// state prolonged (piecewise constant) to the fine grid.
    pub fn fine_state(&self) -> Cow<'_, VoxelState> {
        if self.level == 0 {
            Cow::Borrowed(&self.state)
        } else {
            Cow::Owned(self.prolonged_to(0))
        }
    }

    /// The state prolonged from the current level to the finer `level`.
    fn prolonged_to(&self, level: u32) -> VoxelState {
        let (from, to) = (self.level as usize, level as usize);
        assert!(
            to < from,
            "voxel coarsening: level {to} is not finer than {from}"
        );
        let mut st = prolonged(&self.grids[from - 1], &self.state);
        for g in self.grids[to..from - 1].iter().rev() {
            st = prolonged(g, &st);
        }
        st
    }

    /// Move to `level` with `state` on its grid: the diffusion operator is
    /// rebuilt for that grid and the per-voxel carried steps are forgotten
    /// (the background's is kept).
    fn switch_level(&mut self, level: u32, state: VoxelState) {
        let grid = self.grids[level as usize].clone();
        debug_assert_eq!(state.n_voxels, grid.n_voxels());
        self.diffusion = Diffusion::new(grid, self.diffusion.coeffs().to_vec());
        self.carry.reset_voxels(state.n_voxels);
        self.state = state;
        self.level = level;
        self.level_history.push((self.t, level));
    }

    /// Refine to `level` if the run is coarser.
    fn refine_to(&mut self, level: u32) {
        if self.level > level {
            let st = self.prolonged_to(level);
            self.switch_level(level, st);
        }
    }

    /// Drop one level if the maximum allows it and the state passes
    /// [`can_coarsen`] with the control tolerances.
    fn try_coarsen(&mut self) {
        if self.level >= self.max_coarsen_levels {
            return;
        }
        let grid = &self.grids[self.level as usize];
        let c = &self.control;
        if can_coarsen(grid, &self.state, c.atol, c.rtol, self.coarsen_tol_factor) {
            let st = restricted(grid, &self.state);
            self.switch_level(self.level + 1, st);
        }
    }

    /// Fine voxels per voxel of the current level (`8^level`, the volume ratio).
    fn fine_per_voxel(&self) -> f64 {
        8f64.powi(self.level as i32)
    }

    /// Add `amount` [µmol/L on the fine grid, same layout as the fine state]
    /// and damp the next diffusion step (Rannacher). A coarsened run is first
    /// refined to level 0. dt restarts from `dt_initial` and the carried rsoda
    /// steps are forgotten.
    pub fn inject(&mut self, amount: &VoxelState) {
        self.refine_to(0);
        assert_eq!(amount.data.len(), self.state.data.len());
        for (s, a) in self.state.data.iter_mut().zip(&amount.data) {
            *s += a;
        }
        self.rannacher_pending = true;
        self.dt = self.control.dt_initial;
        self.carry.reset();
        self.record_totals();
    }

    /// Start recording domain sums (records the current state first).
    pub fn track_totals(&mut self) {
        self.totals = Some(Vec::new());
        self.background_totals = self.background.as_ref().map(|_| Vec::new());
        self.record_totals();
    }

    /// Domain sum of every species [µmol/L summed over fine voxels]. At a
    /// coarse level each voxel counts for the `8^level` fine voxels it covers
    /// (its volume), so the sums do not depend on the level.
    pub fn domain_sums(&self) -> Vec<f64> {
        let w = self.fine_per_voxel();
        (0..self.state.n_species)
            .map(|s| self.state.field(s).iter().sum::<f64>() * w)
            .collect()
    }

    fn record_totals(&mut self) {
        if self.totals.is_some() {
            let sums = self.domain_sums();
            if let Some(v) = self.totals.as_mut() {
                v.push((self.t, sums));
            }
            let n = self.state.n_voxels as f64 * self.fine_per_voxel();
            if let (Some(v), Some(bg)) = (self.background_totals.as_mut(), &self.background) {
                v.push(bg.iter().map(|b| b * n).collect());
            }
        }
    }

    fn strang(
        &mut self,
        state: &mut VoxelState,
        bg: &mut Option<Vec<f64>>,
        carry: &mut CarriedSteps,
        t: f64,
        dt: f64,
        damped: bool,
    ) -> Result<()> {
        self.half_reaction(state, bg, carry, t, dt / 2.0)?;
        self.diffuse_shifted(state, bg, dt, damped);
        self.half_reaction(state, bg, carry, t + dt / 2.0, dt / 2.0)
    }

    /// Two consecutive Strang steps of `dt / 2` each, starting at `t`, with the
    /// reaction half-steps at their junction fused into one solve over `dt / 2`
    /// (5 reaction solves instead of 6). `damped` applies to the first diffusion only.
    fn strang_pair(
        &mut self,
        state: &mut VoxelState,
        bg: &mut Option<Vec<f64>>,
        carry: &mut CarriedSteps,
        t: f64,
        dt: f64,
        damped: bool,
    ) -> Result<()> {
        let q = dt / 4.0;
        self.half_reaction(state, bg, carry, t, q)?;
        self.diffuse_shifted(state, bg, dt / 2.0, damped);
        self.half_reaction(state, bg, carry, t + q, 2.0 * q)?;
        self.diffuse_shifted(state, bg, dt / 2.0, false);
        self.half_reaction(state, bg, carry, t + 3.0 * q, q)
    }

    fn diffuse_shifted(
        &mut self,
        state: &mut VoxelState,
        bg: &Option<Vec<f64>>,
        dt: f64,
        damped: bool,
    ) {
        let shift = if self.skip_background {
            bg.as_deref()
        } else {
            None
        };
        self.diffuse(state, shift, dt, damped);
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
        carry: &mut CarriedSteps,
        t: f64,
        h: f64,
    ) -> Result<()> {
        let (mat, o) = (&self.mat, &self.solve);
        let n_bg = if self.skip_background {
            reaction_step_background(mat, state, bg.as_mut(), carry, t, h, o)?
        } else {
            if let Some(b) = bg.as_mut() {
                integrate_one(mat, b, &mut carry.background, t, h, o)?;
            }
            reaction_step_background(mat, state, None, carry, t, h, o)?
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

    /// Advance to exactly `t_target`, checking after each accepted step
    /// whether the grid may drop one coarsening level.
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
            let mut big_carry = self.carry.clone();
            self.strang(&mut big, &mut big_bg, &mut big_carry, self.t, dt, damped)?;
            let mut small = self.state.clone();
            let mut small_bg = self.background.clone();
            let mut small_carry = self.carry.clone();
            let (sb, sc) = (&mut small_bg, &mut small_carry);
            self.strang_pair(&mut small, sb, sc, self.t, dt, damped)?;
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
                self.carry = small_carry;
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
                // Merging only smooths and implicit ADI has no CFL limit: dt is kept.
                self.try_coarsen();
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

    #[cfg(feature = "parallel")]
    #[test]
    fn seeded_run_does_not_depend_on_thread_count() {
        let pool = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
        };
        let one = pool(1).install(|| seeded_run(true));
        let four = pool(4).install(|| seeded_run(true));
        assert_eq!(one.state.data, four.state.data);
        assert_eq!(one.background, four.background);
        assert_eq!(one.steps, four.steps);
        assert_eq!(one.carry, four.carry);
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

    /// `A + B -> C` (k = 1e10 /M/s) with every species dynamic.
    fn ab_env() -> SimEnv {
        let cfg = Config::from_toml_str(
            "[parameters]\n[concentrations.initial]\nA = 0.0\nB = 0.0\nC = 0.0\n[[reactions.k_reaction]]\nreaction = \"A + B -> C\"\nk_value = 1e10\n",
            None,
        )
        .unwrap();
        SimEnv::from_config(&cfg).unwrap()
    }

    /// Indices of A, B and C in the state.
    fn abc(env: &SimEnv) -> [usize; 3] {
        ["A", "B", "C"].map(|s| env.species.index_of(s).unwrap())
    }

    const AB_DIMS: [usize; 3] = [8, 4, 4];

    /// A run of `ab_env` from `state` on an 8x4x4 grid of 5 nm voxels with up to
    /// `max_levels` coarsening levels. `dt_max` keeps steps coming after the grid
    /// has coarsened.
    fn ab_voxel_run(state: VoxelState, max_levels: u32) -> VoxelRun {
        let env = ab_env();
        let [ia, ib, ic] = abc(&env);
        let grid = VoxelGrid::new(AB_DIMS, [5e-9; 3], [0.0; 3]).unwrap();
        let mut coeffs = vec![0.0; env.species.n_dyn()];
        coeffs[ia] = 2e-9;
        coeffs[ib] = 1e-9;
        coeffs[ic] = 1.5e-9;
        let mut run = VoxelRun::new(
            SimMatrices::build(&env).unwrap(),
            Diffusion::new(grid, coeffs),
            state,
            0.0,
            StepControl {
                dt_max: 1e-7,
                ..StepControl::default()
            },
            opts(),
        );
        run.set_coarsening(max_levels, 0.1).unwrap();
        run
    }

    /// A (300 µM) in the low-x half, B (100 µM) in the high-x half: they mix,
    /// react and homogenise. Run to `t_end` with up to `max_levels` levels.
    fn homogenising_run(max_levels: u32, t_end: f64) -> VoxelRun {
        let env = ab_env();
        let [ia, ib, _] = abc(&env);
        let grid = VoxelGrid::new(AB_DIMS, [5e-9; 3], [0.0; 3]).unwrap();
        let mut st = VoxelState::zeros(env.species.n_dyn(), grid.n_voxels());
        for v in 0..grid.n_voxels() {
            if grid.coords(v)[0] < 4 {
                st.field_mut(ia)[v] = 300.0;
            } else {
                st.field_mut(ib)[v] = 100.0;
            }
        }
        let mut run = ab_voxel_run(st, max_levels);
        run.track_totals();
        run.advance_to(t_end).unwrap();
        run
    }

    #[test]
    fn homogenising_run_coarsens_and_keeps_the_totals() {
        let env = ab_env();
        let [ia, ib, ic] = abc(&env);
        let t_end = 3e-6;
        let fine = homogenising_run(0, t_end);
        let coarse = homogenising_run(2, t_end);
        assert_eq!(fine.level(), 0);
        assert!(fine.level_history.is_empty());
        // One level at a time, after accepted steps, with steps left on level 2.
        let levels: Vec<u32> = coarse.level_history.iter().map(|(_, l)| *l).collect();
        assert_eq!(levels, [1, 2]);
        assert!(
            coarse.level_history[1].0 < t_end - 5e-7,
            "{:?}",
            coarse.level_history
        );
        assert_eq!(coarse.level(), 2);
        assert_eq!(coarse.state.n_voxels, 2);
        assert_eq!(coarse.diffusion.grid().dims, [2, 1, 1]);
        assert_eq!(coarse.carry.voxels.len(), 2);
        assert!(coarse.voxel_steps < fine.voxel_steps);
        // Same domain totals as the run without coarsening.
        for (a, b) in coarse.domain_sums().iter().zip(fine.domain_sums()) {
            assert!((a - b).abs() <= 1e-6 * b.abs(), "{a} vs {b}");
        }
        // A + C and B + C hold through the level changes: a coarse voxel counts for the
        // fine voxels it covers.
        let half = (AB_DIMS.iter().product::<usize>() / 2) as f64;
        for (t, sums) in coarse.totals.as_ref().unwrap() {
            for (s, want) in [(ia, 300.0 * half), (ib, 100.0 * half)] {
                let got = sums[s] + sums[ic];
                assert!(
                    (got - want).abs() <= 1e-9 * want,
                    "t {t:e}: {got} vs {want}"
                );
            }
        }
        // The field on the fine grid stays within the control tolerance of the fine run.
        let c = &coarse.control;
        let field = coarse.fine_state();
        assert_eq!(field.n_voxels, fine.state.n_voxels);
        for (a, b) in field.data.iter().zip(&fine.state.data) {
            assert!((a - b).abs() <= c.atol + c.rtol * b.abs(), "{a} vs {b}");
        }
    }

    #[test]
    fn injection_after_coarsening_refines_and_conserves_the_deposit() {
        let env = ab_env();
        let [ia, ib, ic] = abc(&env);
        let grid = VoxelGrid::new(AB_DIMS, [5e-9; 3], [0.0; 3]).unwrap();
        let n_fine = grid.n_voxels();
        let mut y0 = vec![0.0; env.species.n_dyn()];
        y0[ia] = 100.0;
        let mut run = ab_voxel_run(VoxelState::uniform(&y0, n_fine), 2);
        run.background = Some(y0.clone());
        run.track_totals();
        run.advance_to(1e-10).unwrap();
        // A uniform field drops one level per accepted step, up to the maximum, and its
        // voxels stay background on every level.
        assert_eq!(run.level(), 2);
        assert_eq!(run.state.n_voxels, n_fine / 64);
        assert_eq!(run.background_voxel_steps, run.voxel_steps);
        let bg = run.background.clone().unwrap();
        assert_eq!(run.fine_state().data, VoxelState::uniform(&bg, n_fine).data);
        let bg_totals = run.background_totals.as_ref().unwrap().last().unwrap();
        assert_eq!(bg_totals[ia], bg[ia] * n_fine as f64);
        let before = run.domain_sums();
        assert_eq!(before[ia], bg[ia] * n_fine as f64);

        let mut add = VoxelState::zeros(y0.len(), n_fine);
        let (p, q) = (grid.index(5, 2, 1), grid.index(6, 2, 1));
        add.field_mut(ib)[p] = 100.0;
        add.field_mut(ib)[q] = 40.0;
        run.inject(&add);
        assert_eq!(run.level(), 0);
        assert_eq!(run.diffusion.grid(), &grid);
        assert_eq!(run.carry, CarriedSteps::new(n_fine));
        let levels: Vec<u32> = run.level_history.iter().map(|(_, l)| *l).collect();
        assert_eq!(levels, [1, 2, 0]);
        assert_eq!(run.level_history[2].0, 1e-10);
        // Rannacher start as without coarsening.
        assert!(run.rannacher_pending);
        assert_eq!(run.dt(), run.control.dt_initial);
        // The deposit lands on its own fine voxels, on top of the prolonged state.
        let mut want = VoxelState::uniform(&bg, n_fine);
        want.field_mut(ib)[p] += 100.0;
        want.field_mut(ib)[q] += 40.0;
        assert_eq!(run.state, want);
        let after = run.domain_sums();
        assert_eq!(after[ib] - before[ib], 140.0);
        assert_eq!(after[ia], before[ia]);

        // The deposit is conserved as B + C while it reacts and spreads; the fresh
        // deposit keeps the grid fine.
        run.advance_to(1e-8).unwrap();
        let sums = run.domain_sums();
        assert!(sums[ic] > 0.5, "B must react: C = {}", sums[ic]);
        let deposit = sums[ib] + sums[ic];
        assert!((deposit - 140.0).abs() <= 1e-9 * 140.0, "{deposit}");
        let a_total = sums[ia] + sums[ic];
        assert!(
            (a_total - before[ia]).abs() <= 1e-9 * before[ia],
            "{a_total}"
        );
        assert_eq!(run.level(), 0);
    }

    #[test]
    fn set_coarsening_checks_the_grid_and_refines_above_the_maximum() {
        let env = ab_env();
        let n = env.species.n_dyn();
        let make = |dims| {
            let grid = VoxelGrid::new(dims, [5e-9; 3], [0.0; 3]).unwrap();
            let st = VoxelState::zeros(n, grid.n_voxels());
            let diffusion = Diffusion::new(grid, vec![1e-9; n]);
            let mat = SimMatrices::build(&env).unwrap();
            VoxelRun::new(mat, diffusion, st, 0.0, StepControl::default(), opts())
        };
        let mut run = make([8, 4, 6]);
        let err = run.set_coarsening(2, 0.1).unwrap_err().to_string();
        assert!(err.contains("along z (6)"), "{err}");
        for (levels, factor) in [(1, 0.0), (1, -1.0), (1, f64::NAN), (99, 0.1)] {
            assert!(
                run.set_coarsening(levels, factor).is_err(),
                "{levels} {factor}"
            );
        }
        run.set_coarsening(1, 0.5).unwrap();

        let mut run = make(AB_DIMS);
        run.set_coarsening(2, 0.1).unwrap();
        run.advance_to(1e-11).unwrap();
        assert_eq!(run.level(), 2);
        run.set_coarsening(1, 0.1).unwrap();
        assert_eq!(run.level(), 1);
        assert_eq!(run.state.n_voxels, 16);
        assert_eq!(run.diffusion.grid().dims, [4, 2, 2]);
        run.set_coarsening(0, 0.1).unwrap();
        assert_eq!(run.level(), 0);
        assert_eq!(run.state, VoxelState::zeros(n, 128));
        // Off: the run stays fine.
        run.advance_to(2e-11).unwrap();
        assert_eq!(run.level(), 0);
        let levels: Vec<u32> = run.level_history.iter().map(|(_, l)| *l).collect();
        assert_eq!(levels, [1, 2, 1, 0]);
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
