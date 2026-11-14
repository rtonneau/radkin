//! Tau-leaping against the exact SSA.
//!
//! Tau-leaping fires the reactions with the propensities of the start of
//! each leap, so its mean carries a systematic error proportional to
//! `tau_epsilon` on top of the sampling noise. The tests check both parts:
//! with a small epsilon the two ensembles agree within their sampling error,
//! and with the default epsilon within that plus `epsilon` times the peak of
//! the species.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::matrices::SimMatrices;
use radkin::physics::Beam;
use radkin::ssa::{
    EnsembleOptions, SsaConfig, SsaModel, SsaResult, SsaRng, TauLeap, beam_segments, run_ensemble,
    run_trajectory,
};
use radkin::unit_cell::{UnitCell, thinned_times};

const K: f64 = 4.0;

/// No fast equilibrium: leaping pays at large counts.
const RADICALS: &str = r#"
    [concentrations.initial]
    H_r = 1e-6
    OH_r = 3e-6
    e_aq = 2e-6

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

fn ensemble(
    cell: &UnitCell,
    t_span: (f64, f64),
    volume: f64,
    n: usize,
    tau: Option<f64>,
) -> SsaResult {
    let cfg = SsaConfig {
        volume: Some(volume),
        trajectories: n,
        seed: 77,
        tau_leap: tau.is_some(),
        tau_epsilon: tau.unwrap_or(0.03),
        rate_correction: None,
    };
    let options = EnsembleOptions {
        t_eval: Some(thinned_times(t_span.0, t_span.1, 20)),
        ..Default::default()
    };
    run_ensemble(cell, &cfg, t_span, &options).unwrap()
}

/// Check the means and standard deviations of `tau` against `exact`.
/// `epsilon` is the systematic allowance, as a fraction of the species peak.
fn compare(name: &str, exact: &SsaResult, tau: &SsaResult, epsilon: f64) {
    let molecule = exact.molecule();
    let n = exact.trajectories as f64;
    let n_species = exact.mean[0].len();
    let peak: Vec<f64> = (0..n_species)
        .map(|s| exact.mean.iter().map(|row| row[s]).fold(0.0, f64::max))
        .collect();
    let (mut worst, mut worst_std) = (0.0f64, 0.0f64);
    for i in 0..exact.t.len() {
        for s in 0..n_species {
            let (me, mt) = (exact.mean[i][s], tau.mean[i][s]);
            let (se, st) = (exact.std[i][s], tau.std[i][s]);
            // Combined standard error, with the Poisson one as a floor.
            let var = (se * se + st * st) / n;
            let floor = me.max(molecule) * molecule / n;
            let sem = var.max(floor).sqrt();
            let allowed = K * sem + epsilon * peak[s] + molecule;
            assert!(
                (me - mt).abs() <= allowed,
                "{name}: species {s} at t = {:e} s: tau {mt:e} vs exact {me:e} µmol/L \
                 (combined SEM {sem:e})",
                exact.t[i]
            );
            worst = worst.max((me - mt).abs() / sem);
            // Standard deviations: the sampling error of each is about
            // sd / sqrt(2 (n - 1)).
            let sd_allowed = K * se.max(st) / (n - 1.0).sqrt() + epsilon * se + molecule;
            assert!(
                (se - st).abs() <= sd_allowed,
                "{name}: species {s} at t = {:e} s: std {st:e} vs {se:e} µmol/L",
                exact.t[i]
            );
            if se > 10.0 * molecule {
                worst_std = worst_std.max((se - st).abs() / se);
            }
        }
    }
    let events = tau.n_events as f64 / exact.n_events as f64;
    println!(
        "{name}: worst mean difference {worst:.2} combined SEM, worst std difference \
         {:.1} %, events ratio {events:.4}",
        100.0 * worst_std
    );
    assert!((events - 1.0).abs() < 0.02, "{name}: events ratio {events}");
}

#[test]
fn ssa_tau_leaping_matches_exact_at_large_counts() {
    let cell = UnitCell::from_config(&Config::from_toml_str(RADICALS, None).unwrap()).unwrap();
    // 60 000 to 180 000 molecules per species.
    let (t_span, volume, n) = ((0.0, 1e-3), 1e-13, 64);
    let exact = ensemble(&cell, t_span, volume, n, None);
    // A small epsilon: sampling error only.
    let fine = ensemble(&cell, t_span, volume, n, Some(0.001));
    compare("epsilon 0.001", &exact, &fine, 0.0);
    // The default epsilon: its systematic error is allowed for.
    let default = ensemble(&cell, t_span, volume, n, Some(0.03));
    compare("epsilon 0.03", &exact, &default, 0.03);
}

#[test]
fn ssa_tau_leaping_matches_exact_with_a_pulsed_beam() {
    // Pure water: the fast acid-base equilibria make most steps exact.
    let mut cell = UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap();
    let beam = Beam::pulsed_from_dose_per_pulse(3.0, 1e-3, 2e-6, Some(3.0), None, None).unwrap();
    cell.set_beam(beam, "default", Some(2.5)).unwrap();
    let (t_span, volume, n) = ((0.0, 3.5e-3), 1e-14, 48);
    let exact = ensemble(&cell, t_span, volume, n, None);
    let tau = ensemble(&cell, t_span, volume, n, Some(0.03));
    compare("pulsed water", &exact, &tau, 0.03);
}

#[test]
fn ssa_tau_leaping_sources_respect_the_pulses() {
    // Production only: the yield must be dose x G, and nothing may be
    // produced between pulses (a leap never crosses a beam transition).
    let env = SimEnv::from_config(
        &Config::from_toml_str("[reactions.radiolytic]\ne_aq = 2.6", None).unwrap(),
    )
    .unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    let model = SsaModel::new(&mat, 1e-13).unwrap();
    let e = env.species.index_of("e_aq").unwrap();
    let per_gray = model.sources[0].per_gray;
    let beam = Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
    let seg = beam_segments(Some(&beam), (0.0, 0.05)).unwrap();
    // Mid first pulse, just after it, just before the second, the end.
    let t = [1e-6, 2e-6, 9.9e-3, 0.05];
    let doses = [1.0, 2.0, 2.0, 6.0];
    let n = 200;
    let mut sums = [0.0; 4];
    for stream in 0..n {
        let mut rng = SsaRng::new(8, stream);
        let mut stepper = TauLeap::new(&model);
        let traj = run_trajectory(
            &model,
            &seg,
            vec![0; model.n_dyn],
            &t,
            &mut stepper,
            &mut rng,
        )
        .unwrap();
        assert_eq!(traj.rows[1][e], traj.rows[2][e], "beam off: no production");
        assert!(stepper.n_leaps > 0);
        for (sum, row) in sums.iter_mut().zip(&traj.rows) {
            *sum += row[e] as f64;
        }
    }
    for (sum, dose) in sums.iter().zip(doses) {
        let expected = dose * per_gray;
        let mean = sum / n as f64;
        let sem = (expected / n as f64).sqrt();
        assert!(
            (mean - expected).abs() <= K * sem,
            "dose {dose} Gy: mean {mean} vs {expected} (sem {sem})"
        );
    }
}
