//! The SSA ensemble mean against LSODA in the large-volume limit.
//!
//! Criterion, per species and output time, in molecules:
//!
//! ```text
//! |mean - lsoda| <= K * SEM_eff + BIAS
//! ```
//!
//! - `SEM_eff` is the larger of the sample standard error of the mean and
//!   the Poisson one, `sqrt(lsoda / N)`: a species with a mean far below one
//!   molecule often has a sample SEM of 0.
//! - `BIAS` (one molecule) allows for the finite-volume difference between
//!   the mean of the stochastic model and the rate equations.
//!
//! The seed is fixed and the samplers are platform independent, so the test
//! is deterministic. See `docs/ssa.md` for the measured agreement.

use radkin::physics::Beam;
use radkin::result::OdeResult;
use radkin::ssa::{EnsembleOptions, SsaConfig, run_ensemble};
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};

const K: f64 = 4.0;
const BIAS: f64 = 1.0;
/// 10 µm³: 6022 molecules per µmol/L.
const VOLUME: f64 = 1e-14;
const TRAJECTORIES: usize = 64;
const SEED: u64 = 2026;
/// Output times compared, picked among the LSODA points after `T_MIN`.
const N_TIMES: usize = 12;
/// LSODA skips the first 1e-7 s of a beam ON phase (its yield is added at
/// once): compare well after it.
const T_MIN: f64 = 1e-6;

fn pure_water(beam: Beam) -> UnitCell {
    let mut cell = UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap();
    cell.set_beam(beam, "default", Some(2.5)).unwrap();
    cell
}

/// Indices of `N_TIMES` LSODA points spread over `t >= T_MIN`.
fn pick(t: &[f64]) -> Vec<usize> {
    let idx: Vec<usize> = (0..t.len()).filter(|i| t[*i] >= T_MIN).collect();
    let mut out: Vec<usize> = (0..N_TIMES)
        .map(|k| idx[k * (idx.len() - 1) / (N_TIMES - 1)])
        .collect();
    out.dedup();
    out
}

/// Run the ensemble at the picked LSODA times and check every species.
fn compare(name: &str, cell: &UnitCell, lsoda: &OdeResult, y0: Option<Vec<f64>>, no_beam: bool) {
    let t = lsoda.t();
    let picks = pick(t);
    assert!(picks.len() >= N_TIMES / 2, "{name}: too few LSODA points");
    let t_eval: Vec<f64> = picks.iter().map(|i| t[*i]).collect();
    let cfg = SsaConfig {
        volume: Some(VOLUME),
        trajectories: TRAJECTORIES,
        seed: SEED,
        tau_leap: false,
    };
    let options = EnsembleOptions {
        t_eval: Some(t_eval.clone()),
        y0,
        no_beam,
    };
    let ssa = run_ensemble(cell, &cfg, (t[0].min(0.0), *t.last().unwrap()), &options).unwrap();
    assert!(ssa.n_events > 0);
    let molecule = ssa.molecule();
    let sem = ssa.sem();
    let labels = cell.env.species.sorted_dyn_labels();
    let n = TRAJECTORIES as f64;
    let (mut worst, mut checked, mut populated) = (0.0f64, 0, 0);
    for (k, i) in picks.iter().enumerate() {
        for (s, label) in labels.iter().enumerate() {
            let ode = lsoda.y()[*i][s] / molecule;
            let mean = ssa.mean[k][s] / molecule;
            let sem_eff = (sem[k][s] / molecule).max((ode.max(0.0) / n).sqrt());
            let diff = (mean - ode).abs();
            assert!(
                diff <= K * sem_eff + BIAS,
                "{name}: {label} at t = {:e} s: ssa {mean:.3} vs lsoda {ode:.3} molecules \
                 (SEM {sem_eff:.3})",
                t_eval[k]
            );
            if sem_eff > 0.0 {
                worst = worst.max(diff / sem_eff);
            }
            checked += 1;
            // A count large enough for the comparison to be sharp: a relative
            // SEM below 5 %.
            if ode > 100.0 && sem_eff < 0.05 * ode {
                populated += 1;
            }
        }
    }
    println!(
        "{name}: {checked} checks ({populated} on species above 100 molecules), \
         worst |mean - lsoda| = {worst:.2} SEM, {} events per trajectory",
        ssa.n_events / TRAJECTORIES as u64
    );
    assert!(
        populated >= 3 * N_TIMES,
        "{name}: only {populated} sharp checks"
    );
}

#[test]
fn ssa_validation_reaction_only() {
    // Radicals as after a pulse, no beam: the self-reactions OH + OH, H + H
    // and e_aq + e_aq and the cross reactions set the decay.
    let mut cell = pure_water(Beam::constant(1.0, None, None).unwrap());
    let mut y0 = cell.env.initial_values();
    for (label, umol) in [
        ("e_aq", 3.0),
        ("OH_r", 3.0),
        ("H_r", 0.6),
        ("H2O2", 0.7),
        ("H3O_plus", 3.0),
    ] {
        y0[cell.env.species.index_of(label).unwrap()] += umol;
    }
    let options = RunOptions {
        use_beam: false,
        ..Default::default()
    };
    let lsoda = cell.run((0.0, 1e-3), &options, Some(&y0)).unwrap();
    compare("reaction-only", &cell, &lsoda, Some(y0), true);
}

#[test]
fn ssa_validation_continuous_beam() {
    // 1e4 Gy/s for 1 ms: 10 Gy.
    let mut cell = pure_water(Beam::constant(1e4, None, None).unwrap());
    let lsoda = cell.run((0.0, 1e-3), &RunOptions::default(), None).unwrap();
    compare("continuous beam", &cell, &lsoda, None, false);
}

#[test]
fn ssa_validation_pulsed_beam() {
    // 3 pulses of 3 Gy (1.5e6 Gy/s for 2 µs) every millisecond.
    let beam = Beam::pulsed_from_dose_per_pulse(3.0, 1e-3, 2e-6, Some(3.0), None, None).unwrap();
    let mut cell = pure_water(beam);
    let lsoda = cell
        .chunked_run((0.0, 3.5e-3), &ChunkOptions::default())
        .unwrap();
    compare("pulsed beam", &cell, &lsoda, None, false);
}
