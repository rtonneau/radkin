//! Exact SSA on one volume, reaction-only cases with known solutions.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::consts::AVOGADRO;
use radkin::matrices::SimMatrices;
use radkin::ssa::{Exact, Segment, SsaModel, SsaRng, Trajectory, initial_counts, run_trajectory};

const CHEMISTRY: &str = r#"
    [[reactions.k_reaction]]
    reaction = "H_r -> H2"
    k_value = 3.0e3

    [[reactions.k_reaction]]
    reaction = "2 OH_r -> H2O2"
    k_value = 5.5e9
"#;

const VOLUME: f64 = 1e-15;

struct Case {
    env: SimEnv,
    model: SsaModel,
}

fn case() -> Case {
    let env = SimEnv::from_config(&Config::from_toml_str(CHEMISTRY, None).unwrap()).unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    let model = SsaModel::new(&mat, VOLUME).unwrap();
    Case { env, model }
}

impl Case {
    fn idx(&self, label: &str) -> usize {
        self.env.species.index_of(label).unwrap()
    }

    fn counts(&self, pairs: &[(&str, i64)]) -> Vec<i64> {
        let mut n = vec![0; self.model.n_dyn];
        for (l, v) in pairs {
            n[self.idx(l)] = *v;
        }
        n
    }

    fn run(&self, n0: &[i64], t_eval: &[f64], seed: u64, stream: u64) -> Trajectory {
        let seg = [Segment {
            start: 0.0,
            stop: *t_eval.last().unwrap(),
            dose_rate: 0.0,
        }];
        let mut rng = SsaRng::new(seed, stream);
        run_trajectory(&self.model, &seg, n0.to_vec(), t_eval, &mut Exact, &mut rng).unwrap()
    }
}

/// Mean and standard error of the mean of `values`.
fn mean_sem(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    (mean, (var / n).sqrt())
}

#[test]
fn ssa_same_seed_and_stream_give_the_same_trajectory() {
    let c = case();
    let n0 = c.counts(&[("H_r", 500), ("OH_r", 800)]);
    let t = [1e-5, 1e-4, 1e-3];
    let a = c.run(&n0, &t, 42, 3);
    let b = c.run(&n0, &t, 42, 3);
    assert_eq!(a, b);
    assert!(a.n_events > 0);
    assert_ne!(a, c.run(&n0, &t, 42, 4), "another stream");
    assert_ne!(a, c.run(&n0, &t, 43, 3), "another seed");
}

#[test]
fn ssa_first_order_decay_follows_the_binomial_law() {
    // Each of the n0 molecules survives to t with probability exp(-k t).
    let c = case();
    let (h, h2) = (c.idx("H_r"), c.idx("H2"));
    let n0 = c.counts(&[("H_r", 300)]);
    let t = [1e-4, 3e-4, 1e-3];
    let runs: Vec<Trajectory> = (0..400).map(|i| c.run(&n0, &t, 7, i)).collect();
    for (i, ti) in t.iter().enumerate() {
        let left: Vec<f64> = runs.iter().map(|r| r.rows[i][h] as f64).collect();
        for r in &runs {
            assert_eq!(r.rows[i][h] + r.rows[i][h2], 300, "H is conserved");
        }
        let p = (-3.0e3 * ti).exp();
        let (mean, sem) = mean_sem(&left);
        let exact_sem = (300.0 * p * (1.0 - p) / 400.0).sqrt();
        assert!(
            (mean - 300.0 * p).abs() <= 4.0 * exact_sem,
            "t = {ti}: mean {mean} vs {}",
            300.0 * p
        );
        assert!(
            (sem - exact_sem).abs() <= 0.25 * exact_sem,
            "t = {ti}: sem {sem} vs {exact_sem}"
        );
    }
}

#[test]
fn ssa_self_reaction_follows_the_rate_equation_at_large_counts() {
    // Rate law r = k [A]^2 with stoichiometry -2, and k = k_value / 2 for a
    // scalar k_value: d[A]/dt = -k_value [A]^2, n(t) = n0 / (1 + k_value n0 t / (N_A V)).
    let c = case();
    let (oh, h2o2) = (c.idx("OH_r"), c.idx("H2O2"));
    let n_start = 4000;
    let n0 = c.counts(&[("OH_r", n_start)]);
    let kv = 5.5e9 / (AVOGADRO * VOLUME);
    let t = [2e-6, 1e-5, 5e-5];
    let runs: Vec<Trajectory> = (0..200).map(|i| c.run(&n0, &t, 11, i)).collect();
    for (i, ti) in t.iter().enumerate() {
        for r in &runs {
            assert_eq!(r.rows[i][oh] + 2 * r.rows[i][h2o2], n_start, "OH pairs");
        }
        let left: Vec<f64> = runs.iter().map(|r| r.rows[i][oh] as f64).collect();
        let (mean, sem) = mean_sem(&left);
        let ode = n_start as f64 / (1.0 + kv * n_start as f64 * ti);
        // The finite-volume correction to the mean is of the order of one molecule.
        assert!(
            (mean - ode).abs() <= 4.0 * sem + 1.0,
            "t = {ti}: mean {mean} vs {ode} (sem {sem})"
        );
        // A factor 2 error in the propensity would move the mean by hundreds of sem.
        let wrong = n_start as f64 / (1.0 + 0.5 * kv * n_start as f64 * ti);
        assert!((mean - wrong).abs() > 50.0 * sem, "t = {ti}");
    }
}

#[test]
fn ssa_initial_counts_keep_the_mean_concentration() {
    let c = case();
    let k = c.model.molecules_per_umol();
    let mut y = vec![0.0; c.model.n_dyn];
    y[c.idx("H_r")] = 12.3 / k;
    y[c.idx("OH_r")] = 5.0 / k;
    let mut rng = SsaRng::new(5, 0);
    let (h, oh) = (c.idx("H_r"), c.idx("OH_r"));
    let mut sum = 0.0;
    let n = 20_000;
    for _ in 0..n {
        let counts = initial_counts(&c.model, &y, &mut rng);
        assert!(counts[h] == 12 || counts[h] == 13);
        assert!((counts[oh] - 5).abs() <= 1, "5.0 up to round-off");
        sum += counts[h] as f64;
    }
    // Bernoulli(0.3): sd of the mean = sqrt(0.21 / n) = 3.2e-3.
    assert!((sum / n as f64 - 12.3).abs() < 4.0 * 3.3e-3);
}

#[test]
fn ssa_output_times_must_be_sorted_and_inside_the_run() {
    let c = case();
    let n0 = c.counts(&[("H_r", 10)]);
    let seg = [Segment {
        start: 0.0,
        stop: 1.0,
        dose_rate: 0.0,
    }];
    let mut rng = SsaRng::new(1, 0);
    for bad in [vec![0.5, 0.1], vec![2.0], vec![-1.0]] {
        assert!(
            run_trajectory(&c.model, &seg, n0.clone(), &bad, &mut Exact, &mut rng).is_err(),
            "{bad:?}"
        );
    }
}
