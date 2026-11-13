//! Reaction step of the splitting: every voxel is an independent
//! homogeneous cell integrated with rsoda (no beam source; species come from
//! seeds).
//!
//! With the `parallel` feature (default) the voxels are spread over rayon
//! threads. Each voxel's integration is deterministic and independent, so the
//! result is bit-identical to the serial path.

use crate::error::Result;
use crate::matrices::SimMatrices;
use crate::solver::{DoseRateSource, SolveOptions, solve};
use crate::voxel::grid::VoxelState;

/// Integrate every voxel of `state` (dynamic species, µmol/L) over
/// `[t, t + h]`. `opts.t_eval` is ignored.
pub fn reaction_step(
    mat: &SimMatrices,
    state: &mut VoxelState,
    t: f64,
    h: f64,
    opts: &SolveOptions,
) -> Result<()> {
    if h <= 0.0 {
        return Ok(());
    }
    let n_sp = state.n_species;
    let mut aos = vec![0.0; n_sp * state.n_voxels];
    for (v, y) in aos.chunks_mut(n_sp).enumerate() {
        state.gather(v, y);
    }
    let mut o = opts.clone();
    o.t_eval = Some(vec![t + h]);
    if o.first_step >= h {
        o.first_step = h / 1e3;
    }
    let one = |y: &mut [f64]| -> Result<()> {
        let sol = solve(mat, &DoseRateSource::None, (t, t + h), y, &o)?;
        y.copy_from_slice(&sol.y[0]);
        Ok(())
    };
    run_voxels(&mut aos, n_sp, one)?;
    for (v, y) in aos.chunks(n_sp).enumerate() {
        state.scatter(v, y);
    }
    Ok(())
}

#[cfg(feature = "parallel")]
fn run_voxels<F>(aos: &mut [f64], n_sp: usize, f: F) -> Result<()>
where
    F: Fn(&mut [f64]) -> Result<()> + Sync + Send,
{
    use rayon::prelude::*;
    aos.par_chunks_mut(n_sp).try_for_each(f)
}

#[cfg(not(feature = "parallel"))]
fn run_voxels<F>(aos: &mut [f64], n_sp: usize, f: F) -> Result<()>
where
    F: Fn(&mut [f64]) -> Result<()>,
{
    aos.chunks_mut(n_sp).try_for_each(f)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::consts::data_dir;
    use crate::solver::JacobianMode;

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
}
