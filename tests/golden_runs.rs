//! Full runs against radiopyo's golden data (tests/golden/data of radiopyo,
//! exported by tools/export_golden.py). Tolerances are radiopyo's own
//! (`tests/golden/test_golden.py`, `test_paper_tonneau_2025.py`).

use radkin::physics::Beam;
use radkin::result::OdeResult;
use radkin::solver::JacobianMode;
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};
use serde_json::Value;

const RTOL_RUN: f64 = 1e-6;
const ATOL_SPECIES_RUN: f64 = 1e-10;
const ATOL_RATES_RUN_FRAC_OF_MAX: f64 = 1e-7;
const N_POINTS_RELATIVE_SPREAD: f64 = 0.15;

/// Every run is checked with the dense analytic Jacobian (radiopyo's
/// setting) and with the sparse one, analytic and finite-difference.
const MODES: [(JacobianMode, bool, &str); 3] = [
    (JacobianMode::Dense, true, "dense"),
    (JacobianMode::Sparse, true, "sparse"),
    (JacobianMode::Sparse, false, "sparse-fd"),
];

fn mode_options(jacobian: JacobianMode, use_jac: bool) -> RunOptions {
    RunOptions {
        jacobian,
        use_jac,
        ..Default::default()
    }
}

fn runs() -> Value {
    serde_json::from_str(&std::fs::read_to_string("tests/golden/runs.json").unwrap()).unwrap()
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

/// Rows (label -> values at `indices`) of species and rates, labelled as radiopyo's golden.
type Rows = Vec<(String, Vec<f64>)>;

fn result_rows(res: &OdeResult, indices: &[usize]) -> (Rows, Rows) {
    let (labels, rows) = res.concentrations();
    let species = labels
        .iter()
        .enumerate()
        .map(|(j, l)| (l.clone(), indices.iter().map(|&i| rows[i][j]).collect()))
        .collect();
    let (_, rrows) = res.reaction_rates();
    let strs = res.env.reaction_strings();
    let rates = strs
        .iter()
        .enumerate()
        .map(|(j, l)| (l.clone(), indices.iter().map(|&i| rrows[i][j]).collect()))
        .collect();
    (species, rates)
}

fn assert_rows(
    actual: &[(String, Vec<f64>)],
    labels: &[String],
    golden: &Value,
    atol: f64,
    frac: f64,
    what: &str,
) {
    let golden_rows: Vec<Vec<f64>> = golden.as_array().unwrap().iter().map(floats).collect();
    assert_eq!(actual.len(), labels.len(), "{what}: label count");
    let peak = golden_rows
        .iter()
        .flatten()
        .fold(0.0f64, |m, v| m.max(v.abs()));
    let floor = atol + frac * peak;
    let mut failures = Vec::new();
    for (label, exp) in labels.iter().zip(&golden_rows) {
        let got = &actual
            .iter()
            .find(|(l, _)| l == label)
            .unwrap_or_else(|| panic!("{what}: no {label}"))
            .1;
        assert_eq!(got.len(), exp.len(), "{what} {label}");
        for (k, (g, e)) in got.iter().zip(exp).enumerate() {
            if (g - e).abs() > floor + RTOL_RUN * e.abs() {
                failures.push(format!("{label}[{k}]: {g:e} vs {e:e}"));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{what}: {} mismatch(es):\n{}",
        failures.len(),
        failures[..failures.len().min(15)].join("\n")
    );
}

fn assert_run(res: &OdeResult, indices: &[usize], g: &Value, what: &str) {
    let times: Vec<f64> = indices.iter().map(|&i| res.t()[i]).collect();
    for (a, e) in times.iter().zip(floats(&g["time"])) {
        assert!((a - e).abs() <= 1e-12 * e.abs(), "{what} time {a} vs {e}");
    }
    let (species, rates) = result_rows(res, indices);
    assert_rows(
        &species,
        &strings(&g["species_labels"]),
        &g["species"],
        ATOL_SPECIES_RUN,
        0.0,
        &format!("{what} species"),
    );
    assert_rows(
        &rates,
        &strings(&g["rate_labels"]),
        &g["rate"],
        0.0,
        ATOL_RATES_RUN_FRAC_OF_MAX,
        &format!("{what} rates"),
    );
}

fn assert_n_points(res: &OdeResult, g: &Value, what: &str) {
    let n_new: usize = res.n_points.iter().sum();
    let n_ref: f64 = floats(&g["n_points"]).iter().sum();
    assert!(
        (n_new as f64 - n_ref).abs() <= N_POINTS_RELATIVE_SPREAD * n_ref,
        "{what}: {n_ref} -> {n_new} solver points"
    );
}

#[test]
fn labarbe_run_matches_golden() {
    let golden = runs();
    for (jacobian, use_jac, mode) in MODES {
        let mut cell = UnitCell::from_toml("data/config_Labarbe.toml").unwrap();
        let res = cell
            .run((1e-7, 10.0), &mode_options(jacobian, use_jac), None)
            .unwrap();
        let what = format!("labarbe_run/{mode}");
        assert_n_points(&res, &golden["labarbe_run"], &what);
        let last = res.t().len() - 1;
        assert_run(&res, &[0, last], &golden["labarbe_run"], &what);
    }
}

#[test]
fn labarbe_curve_matches_golden() {
    let golden = runs();
    for (jacobian, use_jac, mode) in MODES {
        let mut cell = UnitCell::from_toml("data/config_Labarbe.toml").unwrap();
        let opts = RunOptions {
            t_eval: Some(floats(&golden["labarbe_curve_t_eval"])),
            ..mode_options(jacobian, use_jac)
        };
        let res = cell.run((1e-7, 10.0), &opts, None).unwrap();
        let all: Vec<usize> = (0..res.t().len()).collect();
        assert_run(
            &res,
            &all,
            &golden["labarbe_curve"],
            &format!("labarbe_curve/{mode}"),
        );
    }
}

#[test]
fn chunked_pulsed_run_matches_golden() {
    let golden = runs();
    let g = &golden["chunked_run"];
    for (jacobian, use_jac, mode) in MODES {
        chunked_pulsed_run(g, jacobian, use_jac, mode);
    }
}

fn chunked_pulsed_run(g: &Value, jacobian: JacobianMode, use_jac: bool, mode: &str) {
    let what = format!("chunked_run/{mode}");
    let mut cell = UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap();
    let beam = Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
    cell.set_beam(beam, "default", Some(2.5)).unwrap();
    let opts = ChunkOptions {
        run: mode_options(jacobian, use_jac),
        max_step_size_on: Some(1e-7),
        max_step_size_off: Some(1e-3),
        ..Default::default()
    };
    let res = cell.chunked_run((0.0, 0.05), &opts).unwrap();
    assert_eq!(res.n_points.len(), g["n_points"].as_array().unwrap().len());
    assert_n_points(&res, g, &what);
    let mut edges = std::collections::BTreeSet::new();
    let mut first = 0;
    for n in &res.n_points {
        edges.insert(first);
        edges.insert(first + n - 1);
        first += n;
    }
    let edges: Vec<usize> = edges.into_iter().collect();
    assert_run(&res, &edges, g, &what);
    for (label, ge) in strings(&g["g_labels"]).iter().zip(floats(&g["g_values"])) {
        let r = cell
            .env
            .reactions
            .radiolytic
            .iter()
            .find(|r| &r.label() == label)
            .unwrap();
        assert!((r.ge() - ge).abs() <= 1e-12 * ge, "{what} {label}");
    }
}

#[test]
fn paper_points_match_golden() {
    paper_points(None);
}

/// `data/tonneau_2025.toml` carries the inputs paper_runs sets in code for the
/// "original" scenario (variable G-values, LET 0.56): only the O2 and the beam's
/// dose rate are set here, and the LET comes from the config beam.
#[test]
fn tonneau_2025_config_matches_original_paper_points() {
    let cell = UnitCell::from_toml("data/tonneau_2025.toml").unwrap();
    assert!(cell.env.variable_gvalues);
    let (beam, let_) = cell.beams.get("default").unwrap();
    assert_eq!((beam.dose_rate_peak, *let_), (1e7, Some(0.56)));
    paper_points(Some("data/tonneau_2025.toml"));
}

/// Re-run the paper points of tests/golden/paper_points.json: every point with
/// its own config, or with `config` only the "original" points, without setting
/// the G-value mode or the LET in code.
fn paper_points(config: Option<&str>) {
    let golden: Value =
        serde_json::from_str(&std::fs::read_to_string("tests/golden/paper_points.json").unwrap())
            .unwrap();
    let rtol = golden["rtol"].as_f64().unwrap();
    let mut failures = Vec::new();
    let points: Vec<&Value> = golden["points"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|p| config.is_none() || p["scenario"] == "original")
        .collect();
    assert!(!points.is_empty());
    for (jacobian, use_jac, mode) in MODES {
        for p in &points {
            let mut cell = match config {
                Some(c) => UnitCell::from_toml(c).unwrap(),
                None => {
                    let c = format!("data/{}", p["config"].as_str().unwrap());
                    let mut cell = UnitCell::from_toml(&c).unwrap();
                    cell.env.variable_gvalues = true;
                    cell
                }
            };
            cell.env
                .set_initial_cc("O2", p["o2_initial_um"].as_f64().unwrap() * 1e-6)
                .unwrap();
            let (dr, dose) = (
                p["dose_rate"].as_f64().unwrap(),
                p["dose"].as_f64().unwrap(),
            );
            cell.set_beam(
                Beam::single_pulse(dr, dose, None).unwrap(),
                "default",
                if config.is_some() {
                    None
                } else {
                    p["let"].as_f64()
                },
            )
            .unwrap();
            let opts = ChunkOptions {
                run: mode_options(jacobian, use_jac),
                max_step_size_on: p["max_step_size_on"].as_f64(),
                max_step_size_final: p["max_step_size_final"].as_f64(),
                ..Default::default()
            };
            let res = cell
                .single_pulse_run((0.0, p["t_end"].as_f64().unwrap()), &opts)
                .unwrap();
            let (labels, rows) = res.concentrations();
            let col = |l: &str| labels.iter().position(|x| x == l).unwrap();
            let last = rows.last().unwrap();
            let ph_min = rows
                .iter()
                .map(|r| -(r[col("H3O_plus")] * 1e-6).log10())
                .fold(f64::INFINITY, f64::min);
            for (key, got) in [
                ("o2_final", last[col("O2")]),
                ("h2o2_final", last[col("H2O2")]),
                ("ph_min", ph_min),
            ] {
                let exp = p[key].as_f64().unwrap();
                if (got - exp).abs() > rtol * exp.abs() {
                    failures.push(format!(
                        "{mode} {} dr={dr} pO2={} {key}: {got:e} vs {exp:e} (rel {:.2e})",
                        p["scenario"],
                        p["po2"],
                        (got - exp).abs() / exp.abs()
                    ));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
