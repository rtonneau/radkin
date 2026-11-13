use radkin::RadkinError;
use radkin::physics::{Beam, BeamCollection};
use serde_json::Value;

fn table(v: &Value) -> (String, toml::Table) {
    let mut t = toml::Table::new();
    let mut kind = String::new();
    for (k, x) in v.as_object().unwrap() {
        if k == "type" {
            kind = x.as_str().unwrap().to_string();
        } else {
            t.insert(k.clone(), toml::Value::Float(x.as_f64().unwrap()));
        }
    }
    (kind, t)
}

fn close(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() <= 1e-12 * b.abs()
}

#[test]
fn beams_match_radiopyo_golden() {
    let text = std::fs::read_to_string("tests/golden/beams.json").unwrap();
    let golden: Value = serde_json::from_str(&text).unwrap();
    for (name, case) in golden.as_object().unwrap() {
        let (kind, params) = table(&case["params"]);
        let beam = Beam::from_table(&kind, &params).unwrap();
        let s: Vec<f64> = case["scalars"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        let got = [beam.time_at_max_dose, beam.dose_rate, beam.dose_rate_peak];
        for (g, e) in got.iter().zip(&s) {
            assert!(close(*g, *e), "{name} scalars {got:?} != {s:?}");
        }
        for (i, w) in case["windows"].as_array().unwrap().iter().enumerate() {
            let (start, stop) = (w[0].as_f64().unwrap(), w[1].as_f64().unwrap());
            for (key, rows) in [
                ("matrix", beam.make_matrix(stop, start)),
                ("intervals", beam.make_time_intervals(stop, start)),
            ] {
                let exp = case[key][i].as_array().unwrap();
                assert_eq!(rows.len(), exp.len(), "{name} {key} window {i}: {rows:?}");
                for (r, e) in rows.iter().zip(exp) {
                    assert!(
                        close(r[0], e[0].as_f64().unwrap()) && close(r[1], e[1].as_f64().unwrap()),
                        "{name} {key} window {i}: {rows:?} vs {exp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn short_pulse_is_rejected() {
    let err = Beam::pulsed(1e6, 0.01, 5e-7, None, None).unwrap_err();
    assert!(matches!(err, RadkinError::ShortBeamOn(_)), "{err}");
}

#[test]
fn frequency_sets_the_period() {
    let mut t = toml::Table::new();
    for (k, v) in [
        ("frequency", 90.0),
        ("on_time", 2e-6),
        ("dose_rate_peak", 5.7e5),
        ("max_dose", 10.0),
    ] {
        t.insert(k.into(), toml::Value::Float(v));
    }
    let beam = Beam::from_table("pulsed", &t).unwrap();
    assert_eq!(beam.timer.period, 1.0 / 90.0);
    assert!(beam.at(1e-6) > 0.0 && beam.at(1e-3) == 0.0);
}

#[test]
fn merge_keeps_default_and_labelled_beams() {
    let a: toml::Table = toml::from_str("[[constant]]\ndose_rate = 1.0\nLET = 0.2\n").unwrap();
    let b: toml::Table = toml::from_str(
        "[[pulsed]]\nlabel = \"p\"\ndose_rate_peak = 1e6\nperiod = 0.1\non_time = 2e-6\nn_pulse = 3\n",
    )
    .unwrap();
    let merged = BeamCollection::from_config(&a)
        .unwrap()
        .merge(&BeamCollection::from_config(&b).unwrap())
        .unwrap();
    assert_eq!(merged.labels(), vec!["default", "p"]);
    assert_eq!(merged.get("default").unwrap().1, Some(0.2));
}
