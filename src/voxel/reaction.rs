//! Reaction step of the splitting: every voxel is an independent
//! homogeneous cell integrated with rsoda (no beam source; species come from
//! seeds).
//!
//! Each voxel is integrated straight to the end of the half-step
//! ([`integrate_to`]: no stored output, kernel buffers reused per thread),
//! starting with the step size its previous half-step ended with
//! ([`CarriedSteps`]) instead of the configured first step.
//!
//! With the `parallel` feature (default) the voxels are spread over rayon
//! threads. Each voxel's integration is deterministic and depends only on its
//! own state and carried step, so the result is bit-identical to the serial
//! path.

use std::sync::atomic::{AtomicU64, Ordering};

use crate::error::Result;
use crate::matrices::SimMatrices;
use crate::solver::{DoseRateSource, Integration, SolveOptions, SolveScratch, integrate_to};
use crate::voxel::grid::VoxelState;

/// Absolute tolerance of the background match, as a fraction of the solver
/// atol.
pub const BACKGROUND_ATOL_FACTOR: f64 = 1e-3;
/// Relative tolerance of the background match.
pub const BACKGROUND_RTOL: f64 = 1e-12;

/// Whether voxel state `y` matches the background `bg` within
/// `atol + BACKGROUND_RTOL * |bg|` per species.
pub fn is_background(y: &[f64], bg: &[f64], atol: f64) -> bool {
    y.iter()
        .zip(bg)
        .all(|(a, b)| (a - b).abs() <= atol + BACKGROUND_RTOL * b.abs())
}

/// Counters of the reaction step, summed over all threads since the last
/// [`reset_stats`]. Read-only instrumentation for `examples/voxel_bench.rs`; it
/// does not influence the integration.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReactionStats {
    /// Voxels summed over all reaction half-steps (background voxels included).
    pub voxel_steps: u64,
    /// Background voxels among them (they share one integration).
    pub background_voxel_steps: u64,
    /// rsoda solves actually run (active voxels plus background integrations).
    pub solves: u64,
    /// Right-hand-side evaluations summed over these solves.
    pub nfev: u64,
    /// Jacobian evaluations summed over these solves.
    pub njev: u64,
    /// rsoda steps summed over these solves.
    pub nst: u64,
}

static VOXEL_STEPS: AtomicU64 = AtomicU64::new(0);
static BACKGROUND_VOXEL_STEPS: AtomicU64 = AtomicU64::new(0);
static SOLVES: AtomicU64 = AtomicU64::new(0);
static NFEV: AtomicU64 = AtomicU64::new(0);
static NJEV: AtomicU64 = AtomicU64::new(0);
static NST: AtomicU64 = AtomicU64::new(0);

/// Zero the [`ReactionStats`] counters.
pub fn reset_stats() {
    for c in [
        &VOXEL_STEPS,
        &BACKGROUND_VOXEL_STEPS,
        &SOLVES,
        &NFEV,
        &NJEV,
        &NST,
    ] {
        c.store(0, Ordering::Relaxed);
    }
}

/// The [`ReactionStats`] counters since the last [`reset_stats`].
pub fn stats() -> ReactionStats {
    ReactionStats {
        voxel_steps: VOXEL_STEPS.load(Ordering::Relaxed),
        background_voxel_steps: BACKGROUND_VOXEL_STEPS.load(Ordering::Relaxed),
        solves: SOLVES.load(Ordering::Relaxed),
        nfev: NFEV.load(Ordering::Relaxed),
        njev: NJEV.load(Ordering::Relaxed),
        nst: NST.load(Ordering::Relaxed),
    }
}

fn count_solve(r: &Integration) {
    SOLVES.fetch_add(1, Ordering::Relaxed);
    NFEV.fetch_add(r.nfev as u64, Ordering::Relaxed);
    NJEV.fetch_add(r.njev as u64, Ordering::Relaxed);
    NST.fetch_add(r.nst as u64, Ordering::Relaxed);
}

/// Step sizes [s] carried from one reaction half-step to the next: one per
/// voxel and one for the background integration. `0` means none: the
/// half-step starts with the configured first step.
///
/// A hint only (each solve still controls its own error), kept in memory and
/// not stored in checkpoints.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CarriedSteps {
    pub voxels: Vec<f64>,
    pub background: f64,
}

impl CarriedSteps {
    /// No carried step for `n_voxels` voxels and the background.
    pub fn new(n_voxels: usize) -> Self {
        Self {
            voxels: vec![0.0; n_voxels],
            background: 0.0,
        }
    }

    /// Forget every carried step.
    pub fn reset(&mut self) {
        self.voxels.fill(0.0);
        self.background = 0.0;
    }
}

/// Fraction of the carried step taken as the first step of the next
/// half-step. rsoda restarts at order 1, where the step the previous solve
/// ended with (at a higher order) fails the error test; 0.01 gave the fewest
/// right-hand-side calls on the voxel benchmark (`docs/perf/voxel.md`).
pub const CARRY_FACTOR: f64 = 0.01;

/// First step of a half-step of length `h`: `CARRY_FACTOR` times the carried
/// step if any, else `opts.first_step`; `h / 1e3` if that does not fit in `h`.
fn first_step(carried: f64, opts: &SolveOptions, h: f64) -> f64 {
    let f = if carried > 0.0 {
        CARRY_FACTOR * carried
    } else {
        opts.first_step
    };
    if f >= h { h / 1e3 } else { f }
}

/// Integrate `y` over `[t, t + h]` from the step `carried`, then set `carried`
/// to the step it ended with.
fn integrate_carried(
    mat: &SimMatrices,
    y: &mut [f64],
    carried: &mut f64,
    t: f64,
    h: f64,
    opts: &SolveOptions,
    scratch: &mut SolveScratch,
) -> Result<Integration> {
    let r = integrate_to(
        mat,
        &DoseRateSource::None,
        (t, t + h),
        y,
        opts,
        first_step(*carried, opts, h),
        scratch,
    )?;
    *carried = r.h_last;
    Ok(r)
}

/// Integrate every voxel of `state` (dynamic species, µmol/L) over
/// `[t, t + h]`, without carried steps. `opts.t_eval` is ignored.
pub fn reaction_step(
    mat: &SimMatrices,
    state: &mut VoxelState,
    t: f64,
    h: f64,
    opts: &SolveOptions,
) -> Result<()> {
    let mut carry = CarriedSteps::new(state.n_voxels);
    reaction_step_background(mat, state, None, &mut carry, t, h, opts).map(|_| ())
}

/// Integrate one homogeneous state `y` [µmol/L] over `[t, t + h]` (no
/// source), with the options of [`reaction_step`], from the step `carried`
/// (updated, see [`CarriedSteps`]).
pub fn integrate_one(
    mat: &SimMatrices,
    y: &mut [f64],
    carried: &mut f64,
    t: f64,
    h: f64,
    opts: &SolveOptions,
) -> Result<()> {
    if h <= 0.0 {
        return Ok(());
    }
    let mut scratch = SolveScratch::new(mat, opts.jacobian);
    integrate_carried(mat, y, carried, t, h, opts, &mut scratch).map(|_| ())
}

/// [`reaction_step`] with carried steps and background voxels.
///
/// Every integrated voxel starts from its step in `carry.voxels` and leaves
/// the step it ended with there (`carry` is reset if its length is not the
/// number of voxels). When `background` is given, it is integrated once over
/// `[t, t + h]` (from `carry.background`) and every voxel that matches it at
/// the start (see [`is_background`], atol `BACKGROUND_ATOL_FACTOR *
/// opts.atol`) receives that result and that step instead of its own rsoda
/// run. The other voxels are integrated exactly as without a background.
/// Returns the number of background voxels.
pub fn reaction_step_background(
    mat: &SimMatrices,
    state: &mut VoxelState,
    background: Option<&mut Vec<f64>>,
    carry: &mut CarriedSteps,
    t: f64,
    h: f64,
    opts: &SolveOptions,
) -> Result<usize> {
    if h <= 0.0 {
        return Ok(0);
    }
    let n_sp = state.n_species;
    if carry.voxels.len() != state.n_voxels {
        *carry = CarriedSteps::new(state.n_voxels);
    }
    VOXEL_STEPS.fetch_add(state.n_voxels as u64, Ordering::Relaxed);
    let one = |scratch: &mut SolveScratch, y: &mut [f64], carried: &mut f64| -> Result<()> {
        let r = integrate_carried(mat, y, carried, t, h, opts, scratch)?;
        count_solve(&r);
        Ok(())
    };
    let new_scratch = || SolveScratch::new(mat, opts.jacobian);
    let mut y = vec![0.0; n_sp];
    let Some(bg) = background else {
        let mut aos = vec![0.0; n_sp * state.n_voxels];
        for (v, y) in aos.chunks_mut(n_sp).enumerate() {
            state.gather(v, y);
        }
        run_voxels(&mut aos, &mut carry.voxels, n_sp, new_scratch, one)?;
        for (v, y) in aos.chunks(n_sp).enumerate() {
            state.scatter(v, y);
        }
        return Ok(0);
    };
    assert_eq!(bg.len(), n_sp);
    let atol = BACKGROUND_ATOL_FACTOR * opts.atol;
    let mut active = Vec::new();
    let mut aos = Vec::new();
    let mut steps = Vec::new();
    for (v, &step) in carry.voxels.iter().enumerate().take(state.n_voxels) {
        state.gather(v, &mut y);
        if !is_background(&y, bg, atol) {
            active.push(v);
            aos.extend_from_slice(&y);
            steps.push(step);
        }
    }
    let n_background = state.n_voxels - active.len();
    BACKGROUND_VOXEL_STEPS.fetch_add(n_background as u64, Ordering::Relaxed);
    run_voxels(&mut aos, &mut steps, n_sp, new_scratch, one)?;
    // Always advanced, so it stays the background of later steps.
    one(&mut new_scratch(), bg, &mut carry.background)?;
    let mut next = active.iter().zip(aos.chunks(n_sp).zip(&steps)).peekable();
    for v in 0..state.n_voxels {
        match next.peek() {
            Some((a, (y, s))) if **a == v => {
                state.scatter(v, y);
                carry.voxels[v] = **s;
                next.next();
            }
            _ => {
                state.scatter(v, bg);
                carry.voxels[v] = carry.background;
            }
        }
    }
    Ok(n_background)
}

/// Run `f` on every voxel of `aos` (with its carried step) and a scratch of
/// the thread running it, built by `init` on the thread's first voxel.
///
/// One scratch per pool thread, not rayon's `try_for_each_init`: that builds
/// one per split job, which on the benchmark was almost one per voxel.
#[cfg(feature = "parallel")]
fn run_voxels<I, F>(aos: &mut [f64], steps: &mut [f64], n_sp: usize, init: I, f: F) -> Result<()>
where
    I: Fn() -> SolveScratch + Sync + Send,
    F: Fn(&mut SolveScratch, &mut [f64], &mut f64) -> Result<()> + Sync + Send,
{
    use rayon::prelude::*;
    use std::sync::Mutex;
    // Slot `i` is only locked by pool thread `i`, so never contended.
    let scratches: Vec<Mutex<Option<SolveScratch>>> = (0..rayon::current_num_threads())
        .map(|_| Mutex::new(None))
        .collect();
    aos.par_chunks_mut(n_sp)
        .zip(steps.par_iter_mut())
        .try_for_each(|(y, s)| {
            let i = rayon::current_thread_index().unwrap_or(0) % scratches.len();
            let mut slot = scratches[i].lock().expect("scratch slot poisoned");
            f(slot.get_or_insert_with(&init), y, s)
        })
}

#[cfg(not(feature = "parallel"))]
fn run_voxels<I, F>(aos: &mut [f64], steps: &mut [f64], n_sp: usize, init: I, f: F) -> Result<()>
where
    I: Fn() -> SolveScratch,
    F: Fn(&mut SolveScratch, &mut [f64], &mut f64) -> Result<()>,
{
    let mut scratch = init();
    aos.chunks_mut(n_sp)
        .zip(steps.iter_mut())
        .try_for_each(|(y, s)| f(&mut scratch, y, s))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::consts::data_dir;
    use crate::solver::{JacobianMode, solve};

    /// Water radiolysis chemistry with a radiolytic initial burst.
    pub(crate) fn water_env() -> SimEnv {
        let cfg = Config::load(data_dir().join("CONFIG_PAPER_01.toml")).unwrap();
        let mut env = SimEnv::from_config(&cfg).unwrap();
        for (s, c) in [
            ("e_aq", 2e-6),
            ("OH_r", 2e-6),
            ("H3O_plus", 2e-6),
            ("H_r", 4e-7),
            ("H2O2", 5e-7),
        ] {
            env.set_initial_cc(s, c).unwrap();
        }
        env
    }

    pub(crate) fn opts() -> SolveOptions {
        SolveOptions {
            rtol: 1e-10,
            atol: 1e-12,
            first_step: 1e-15,
            max_step: None,
            use_jac: true,
            jacobian: JacobianMode::Dense,
            t_eval: None,
        }
    }

    #[test]
    fn uniform_grid_matches_unit_cell() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let t_end = 1e-3;
        let o = SolveOptions {
            t_eval: Some(vec![t_end]),
            ..opts()
        };
        let reference = solve(&mat, &DoseRateSource::None, (0.0, t_end), &y0, &o).unwrap();
        let reference = &reference.y[0];
        let ie = env.species.index_of("e_aq").unwrap();
        assert!(reference[ie] < 0.5 * y0[ie], "chemistry must actually run");

        let mut st = VoxelState::uniform(&y0, 8);
        // Log-spaced reaction steps.
        let mut t = 0.0;
        for k in 0..12 {
            let next = t_end * 10f64.powi(k - 11);
            reaction_step(&mat, &mut st, t, next - t, &opts()).unwrap();
            t = next;
        }
        let mut y = vec![0.0; y0.len()];
        for v in 0..8 {
            st.gather(v, &mut y);
            for (a, b) in y.iter().zip(reference) {
                let tol = 1e-6 * b.abs().max(1e-6);
                assert!((a - b).abs() <= tol, "voxel {v}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn heterogeneous_voxels_are_independent() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let n = 16;
        let mut st = VoxelState::uniform(&y0, n);
        for v in 0..n {
            let mut y = y0.clone();
            for x in &mut y {
                *x *= 1.0 + v as f64 / 4.0;
            }
            st.scatter(v, &y);
        }
        let mut serial = st.clone();
        reaction_step(&mat, &mut st, 0.0, 1e-6, &opts()).unwrap();
        // Serial reference, one voxel at a time.
        let mut y = vec![0.0; y0.len()];
        for v in 0..n {
            serial.gather(v, &mut y);
            let mut one = VoxelState::uniform(&y, 1);
            reaction_step(&mat, &mut one, 0.0, 1e-6, &opts()).unwrap();
            let mut out = vec![0.0; y.len()];
            one.gather(0, &mut out);
            serial.scatter(v, &out);
        }
        assert_eq!(st.data, serial.data, "parallel and serial differ");
    }

    /// Heterogeneous voxels and the background over three half-steps with
    /// carried steps.
    fn carried_run() -> (VoxelState, Vec<f64>, CarriedSteps) {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        let n = 24;
        let mut st = VoxelState::uniform(&y0, n);
        for v in (0..n).filter(|v| v % 3 != 0) {
            let mut y = y0.clone();
            for x in &mut y {
                *x *= 1.0 + v as f64 / 4.0;
            }
            st.scatter(v, &y);
        }
        let mut bg = y0;
        let mut carry = CarriedSteps::new(n);
        let mut t = 0.0;
        for h in [1e-9, 3e-9, 1e-8] {
            reaction_step_background(&mat, &mut st, Some(&mut bg), &mut carry, t, h, &opts())
                .unwrap();
            t += h;
        }
        (st, bg, carry)
    }

    #[test]
    fn carried_steps_follow_each_voxel() {
        let env = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let (st, bg, carry) = carried_run();
        assert!(carry.voxels.iter().all(|h| *h > 0.0));
        // Serial reference: each voxel alone, with its own carried step.
        let y0 = env.initial_values();
        let mut y = vec![0.0; y0.len()];
        for v in 0..st.n_voxels {
            y.copy_from_slice(&y0);
            if v % 3 != 0 {
                for x in &mut y {
                    *x *= 1.0 + v as f64 / 4.0;
                }
            }
            let mut h_carried = 0.0;
            let mut t = 0.0;
            for h in [1e-9, 3e-9, 1e-8] {
                integrate_one(&mat, &mut y, &mut h_carried, t, h, &opts()).unwrap();
                t += h;
            }
            let mut got = vec![0.0; y.len()];
            st.gather(v, &mut got);
            assert_eq!(got, y, "voxel {v}");
            assert_eq!(carry.voxels[v], h_carried, "voxel {v}");
            if v % 3 == 0 {
                assert_eq!(got, bg, "background voxel {v}");
                assert_eq!(carry.voxels[v], carry.background);
            }
        }
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn carried_steps_do_not_depend_on_thread_count() {
        let pool = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
        };
        let one = pool(1).install(carried_run);
        let four = pool(4).install(carried_run);
        assert_eq!(one, four);
    }
}
