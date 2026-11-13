//! `dense::lsoda_dense_output` against scipy's LSODA dense output on real
//! steps (`tests/golden/dense_output.json`, from `tools/export_golden.py`).

use radkin::dense::lsoda_dense_output;
use rsoda::Nordsieck;
use serde_json::Value;

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn usize_of(v: &Value) -> usize {
    v.as_u64().unwrap() as usize
}

#[test]
fn dense_output_matches_scipy() {
    let golden: Value =
        serde_json::from_str(&std::fs::read_to_string("tests/golden/dense_output.json").unwrap())
            .unwrap();
    let steps = golden["steps"].as_array().unwrap();
    assert!(steps.len() > 100);
    let (mut lowered, mut worst) = (0, 0.0_f64);
    for (k, s) in steps.iter().enumerate() {
        let yh = floats(&s["yh"]);
        let nord = Nordsieck {
            yh: &yh,
            nyh: usize_of(&s["nyh"]),
            n: usize_of(&s["n"]),
            nq: usize_of(&s["nq"]),
            nqu: usize_of(&s["nqu"]),
            tn: s["tn"].as_f64().unwrap(),
            h: s["h"].as_f64().unwrap(),
            hu: s["hu"].as_f64().unwrap(),
        };
        lowered += usize::from(nord.nq < nord.nqu);
        let t_step = s["t_step"].as_f64().unwrap();
        let mut out = vec![0.0; nord.n];
        for (t, want) in floats(&s["t"]).into_iter().zip(s["y"].as_array().unwrap()) {
            let want = floats(want);
            lsoda_dense_output(&nord, t_step, t, &mut out);
            let scale = want.iter().fold(0.0_f64, |m, w| m.max(w.abs()));
            for (i, (&a, &e)) in out.iter().zip(&want).enumerate() {
                let err = (a - e).abs();
                worst = worst.max(err / (e.abs() + 1e-300));
                assert!(
                    err <= 1e-12 * e.abs() + 1e-15 * scale,
                    "step {k}, t = {t:e}, component {i}: {a:e} vs scipy {e:e}"
                );
            }
        }
    }
    assert!(lowered > 0, "the golden data has no order decrease");
    eprintln!("worst relative difference: {worst:e}");
}
