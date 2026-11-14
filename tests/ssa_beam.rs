//! Radiolytic sources of the exact SSA: G-value x dose rate as zero-order
//! source reactions.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::matrices::SimMatrices;
use radkin::physics::Beam;
use radkin::ssa::{
    Exact, SsaModel, SsaRng, Trajectory, beam_segments, delivered_dose, run_trajectory,
};

/// Production only: no reaction consumes the species.
const SOURCES: &str = r#"
    [reactions.radiolytic]
    e_aq = 2.6
    OH_r = 2.7
"#;

const VOLUME: f64 = 1e-15;

struct Case {
    env: SimEnv,
    model: SsaModel,
}

fn case() -> Case {
    let env = SimEnv::from_config(&Config::from_toml_str(SOURCES, None).unwrap()).unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    let model = SsaModel::new(&mat, VOLUME).unwrap();
    Case { env, model }
}

impl Case {
    /// Molecules of `label` produced per gray: yield [mol/L/Gy] * N_A V.
    fn per_gray(&self, label: &str) -> f64 {
        let r = self
            .env
            .reactions
            .radiolytic
            .iter()
            .find(|r| r.species == label)
            .unwrap();
        r.kr * self.model.nav
    }

    fn run(&self, beam: Option<&Beam>, t_eval: &[f64], stream: u64) -> Trajectory {
        let seg = beam_segments(beam, (0.0, *t_eval.last().unwrap())).unwrap();
        let mut rng = SsaRng::new(2026, stream);
        let n0 = vec![0; self.model.n_dyn];
        run_trajectory(&self.model, &seg, n0, t_eval, &mut Exact, &mut rng).unwrap()
    }
}

/// The ensemble mean of `label` at `t_eval[i]` against `dose * yield * N_A V`
/// (a Poisson count: variance = mean).
fn check_yield(c: &Case, runs: &[Trajectory], i: usize, label: &str, dose: f64) {
    let s = c.env.species.index_of(label).unwrap();
    let n = runs.len() as f64;
    let mean = runs.iter().map(|r| r.rows[i][s] as f64).sum::<f64>() / n;
    let var = runs
        .iter()
        .map(|r| (r.rows[i][s] as f64 - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let expected = dose * c.per_gray(label);
    let sem = (expected / n).sqrt();
    assert!(
        (mean - expected).abs() <= 4.0 * sem,
        "{label} at point {i}: mean {mean} vs {expected} (sem {sem})"
    );
    // Sample variance of a Poisson count: relative sd about sqrt(2 / n).
    assert!(
        (var - expected).abs() <= 5.0 * expected * (2.0 / n).sqrt(),
        "{label} at point {i}: variance {var} vs {expected}"
    );
}

#[test]
fn ssa_sources_are_built_from_the_g_values() {
    let c = case();
    assert_eq!(c.model.sources.len(), 2);
    // 2.6 molecules / 100 eV is 2.6 / (100 e) per joule, and 1 Gy in 1e-15 L
    // of water (1e-15 kg) is 1e-15 J.
    let expected = 2.6 / (100.0 * 1.602176634e-19) * 1e-15;
    assert!((c.per_gray("e_aq") - expected).abs() <= 1e-9 * expected);
    let e = c.env.species.index_of("e_aq").unwrap();
    let source = c.model.sources.iter().find(|s| s.species == e).unwrap();
    assert!((source.per_gray - expected).abs() <= 1e-9 * expected);
}

#[test]
fn ssa_continuous_beam_produces_yield_times_dose() {
    let c = case();
    // 50 Gy/s: 0.5, 2 and 5 Gy at the output times.
    let beam = Beam::constant(50.0, None, None).unwrap();
    let t = [0.01, 0.04, 0.1];
    let runs: Vec<Trajectory> = (0..300).map(|i| c.run(Some(&beam), &t, i)).collect();
    for (i, ti) in t.iter().enumerate() {
        for label in ["e_aq", "OH_r"] {
            check_yield(&c, &runs, i, label, 50.0 * ti);
        }
    }
}

#[test]
fn ssa_pulsed_beam_produces_yield_times_dose() {
    let c = case();
    // 3 pulses of 2 Gy (1e6 Gy/s for 2 µs) every 10 ms.
    let beam = Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
    // Mid first pulse, between pulses, after the second pulse, at the end.
    let t = [1e-6, 5e-3, 0.015, 0.05];
    let doses = [1.0, 2.0, 4.0, 6.0];
    let seg = beam_segments(Some(&beam), (0.0, 0.05)).unwrap();
    assert!((delivered_dose(&seg) - 6.0).abs() < 1e-6);
    let runs: Vec<Trajectory> = (0..300).map(|i| c.run(Some(&beam), &t, i)).collect();
    for (i, dose) in doses.iter().enumerate() {
        for label in ["e_aq", "OH_r"] {
            check_yield(&c, &runs, i, label, *dose);
        }
    }
    // Nothing is produced between the pulses.
    let e = c.env.species.index_of("e_aq").unwrap();
    let t2 = [2e-6, 9.9e-3];
    for i in 0..20 {
        let r = c.run(Some(&beam), &t2, i);
        assert_eq!(r.rows[0][e], r.rows[1][e]);
    }
}

#[test]
fn ssa_run_without_beam_has_no_source_event() {
    let c = case();
    let r = c.run(None, &[1e-3, 1.0], 0);
    assert_eq!(r.n_events, 0);
    assert!(r.rows.iter().all(|row| row.iter().all(|n| *n == 0)));
}
