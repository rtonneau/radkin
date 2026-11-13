use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::gvalues::GValuesBoscolo;
use radkin::kernel::{Workspace, derive, jacobian};
use radkin::matrices::SimMatrices;
use radkin::physics::{Beam, BeamCollection};
use serde_json::Value;

struct Case {
    env: SimEnv,
    mat: SimMatrices,
    beam: Beam,
    labels: Vec<String>,
}

fn build(name: &str) -> Case {
    let path = if name == "labarbe" {
        "data/config_Labarbe.toml"
    } else {
        "data/advanced_water_radiolysis.toml"
    };
    let cfg = Config::load(path).unwrap();
    let mut env = SimEnv::from_config(&cfg).unwrap();
    let beams = BeamCollection::from_config(&cfg.beam).unwrap();
    let (beam, let_) = match beams.current() {
        Ok(b) => (b.clone(), beams.current_let().unwrap()),
        Err(_) => (Beam::constant(1e7, Some(10.0), None).unwrap(), 0.201),
    };
    let gv = GValuesBoscolo::load(&radkin::consts::data_dir()).unwrap();
    env.apply_let(&gv, let_, Some(beam.dose_rate_peak)).unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    let labels = env.species.sorted_dyn_labels();
    Case {
        env,
        mat,
        beam,
        labels,
    }
}

fn golden() -> Value {
    serde_json::from_str(&std::fs::read_to_string("tests/golden/derive.json").unwrap()).unwrap()
}

fn vector(state: &Value, labels: &[String]) -> Vec<f64> {
    labels.iter().map(|l| state[l].as_f64().unwrap()).collect()
}

fn assert_close(actual: &[f64], expected: &[f64], what: &str) {
    let peak = expected.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - e).abs() <= 1e-12 * e.abs() + 1e-12 * peak,
            "{what}[{i}]: {a} vs {e}"
        );
    }
}

fn run_derive(c: &Case, y: &[f64], dose_rate: f64) -> Vec<f64> {
    let mut ws = Workspace::new(&c.mat);
    let mut dydt = vec![0.0; y.len()];
    derive(y, dose_rate, &c.mat, &mut ws, &mut dydt);
    dydt
}

fn run_jac(c: &Case, y: &[f64]) -> Vec<f64> {
    let n = y.len();
    let mut ws = Workspace::new(&c.mat);
    let mut pd = vec![0.0; n * n];
    jacobian(y, &c.mat, &mut ws, &mut pd);
    pd
}

/// Golden rows `[i][j]` -> column-major `pd[i + j*n]`.
fn col_major(rows: &Value, n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for (i, row) in rows.as_array().unwrap().iter().enumerate() {
        for (j, v) in row.as_array().unwrap().iter().enumerate() {
            out[i + j * n] = v.as_f64().unwrap();
        }
    }
    out
}

#[test]
fn derive_matches_golden() {
    let g = golden();
    for name in ["labarbe", "pure_water"] {
        let c = build(name);
        assert_eq!(c.labels.len(), c.env.species.n_dyn());
        for (key, expected) in g[name]["dydt"].as_object().unwrap() {
            let (state, case) = key.split_once("__").unwrap();
            let y = vector(&g[name]["states"][state], &c.labels);
            let dr = match case {
                "beam_on" => c.beam.at(5e-7),
                "beam_off" => c.beam.at(1e-3),
                _ => 0.0,
            };
            assert_close(
                &run_derive(&c, &y, dr),
                &vector(expected, &c.labels),
                &format!("{name}/{key}"),
            );
        }
    }
}

#[test]
fn jacobian_matches_golden() {
    let g = golden();
    for name in ["labarbe", "pure_water"] {
        let c = build(name);
        let n = c.labels.len();
        let y = vector(&g[name]["states"]["synthetic"], &c.labels);
        let exp = col_major(&g[name]["jacobian"]["synthetic__beam_on"], n);
        assert_close(&run_jac(&c, &y), &exp, &format!("{name} jacobian"));
    }
}

#[test]
fn overshoot_path_matches_golden() {
    let g = golden();
    for name in ["labarbe", "pure_water"] {
        let c = build(name);
        let n = c.labels.len();
        let y = vector(&g[name]["overshoot_state"], &c.labels);
        let dr = c.beam.at(5e-7);
        assert_close(
            &run_derive(&c, &y, dr),
            &vector(&g[name]["overshoot_dydt"], &c.labels),
            &format!("{name} overshoot dydt"),
        );
        let exp = col_major(&g[name]["overshoot_jacobian"], n);
        assert_close(
            &run_jac(&c, &y),
            &exp,
            &format!("{name} overshoot jacobian"),
        );
    }
}

#[test]
fn jacobian_matches_finite_differences() {
    let g = golden();
    let c = build("labarbe");
    let n = c.labels.len();
    for y in [
        vector(&g["labarbe"]["states"]["synthetic"], &c.labels),
        vector(&g["labarbe"]["overshoot_state"], &c.labels),
    ] {
        let jac = run_jac(&c, &y);
        let peak = jac.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for j in 0..n {
            let h = 1e-4 * y[j].abs().max(1e-3);
            if y[j] < 0.0 && y[j] + h > 0.0 || y[j] >= 0.0 && y[j] - h < 0.0 {
                continue; // kink at 0
            }
            let (mut yp, mut ym) = (y.clone(), y.clone());
            yp[j] += h;
            ym[j] -= h;
            let (fp, fm) = (run_derive(&c, &yp, 0.0), run_derive(&c, &ym, 0.0));
            for i in 0..n {
                let fd = (fp[i] - fm[i]) / (2.0 * h);
                let a = jac[i + j * n];
                assert!(
                    (fd - a).abs() <= 1e-3 * a.abs() + 1e-7 * peak,
                    "J[{i},{j}] fd={fd} analytic={a}"
                );
            }
        }
    }
}
