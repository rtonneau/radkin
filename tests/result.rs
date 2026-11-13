use radkin::unit_cell::{RunOptions, UnitCell};
use serde_json::Value;

fn close(a: f64, b: f64, rtol: f64, atol: f64) -> bool {
    (a - b).abs() <= rtol * b.abs() + atol
}

fn as_map(v: &Value) -> Vec<(String, f64)> {
    v.as_object()
        .unwrap()
        .iter()
        .map(|(k, x)| (k.clone(), x.as_f64().unwrap()))
        .collect()
}

fn check(what: &str, got: &[(String, f64)], exp: &[(String, f64)]) {
    let peak = exp.iter().fold(0.0f64, |m, (_, v)| m.max(v.abs()));
    assert_eq!(got.len(), exp.len(), "{what}: {got:?} vs {exp:?}");
    for (label, e) in exp {
        let g = got
            .iter()
            .find(|(l, _)| l == label)
            .unwrap_or_else(|| panic!("{what}: missing {label}"))
            .1;
        assert!(
            close(g, *e, 1e-6, 1e-9 * peak),
            "{what} {label}: {g} vs {e}"
        );
    }
}

#[test]
fn labarbe_short_run_integrals_match_radiopyo() {
    let golden: Value =
        serde_json::from_str(&std::fs::read_to_string("tests/golden/runs.json").unwrap()).unwrap();
    let g = &golden["labarbe_short_integrals"];
    let t_eval: Vec<f64> = g["t_eval"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    let mut cell = UnitCell::from_toml("data/config_Labarbe.toml").unwrap();
    let res = cell
        .run(
            (1e-7, 1e-3),
            &RunOptions {
                t_eval: Some(t_eval),
                ..Default::default()
            },
            None,
        )
        .unwrap();
    let (labels, rows) = res.concentrations();
    let golden_species = g["species"].as_object().unwrap();
    assert_eq!(labels, golden_species.keys().cloned().collect::<Vec<_>>());
    for (j, label) in labels.iter().enumerate() {
        let exp: Vec<f64> = golden_species[label]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        for (i, row) in rows.iter().enumerate() {
            assert!(
                close(row[j], exp[i], 1e-6, 1e-10),
                "{label}[{i}]: {} vs {}",
                row[j],
                exp[i]
            );
        }
    }
    check(
        "integrate_species",
        &res.integrate_species(0.0, f64::INFINITY),
        &as_map(&g["integrate_species"]),
    );
    check(
        "production H2O2",
        &res.integrate_production_of("H2O2", None, None, true, true)
            .unwrap(),
        &as_map(&g["production_H2O2"]),
    );
    check(
        "consumption O2",
        &res.integrate_consumption_of("O2", None, None, true)
            .unwrap(),
        &as_map(&g["consumption_O2"]),
    );
    let prod = res
        .integrate_production_of("H2O2", None, None, true, true)
        .unwrap();
    assert!(
        prod.windows(2).all(|w| w[0].1 >= w[1].1),
        "sorted descending"
    );
}

#[test]
fn rate_columns_follow_radiopyo_order() {
    let mut cell = UnitCell::from_toml("data/config_Labarbe.toml").unwrap();
    let res = cell
        .run((1e-7, 1e-5), &RunOptions::default(), None)
        .unwrap();
    let (labels, rows) = res.reaction_rates();
    assert_eq!(labels, cell.env.reaction_labels());
    assert_eq!(rows.len(), res.t().len());
    assert!(labels.last().unwrap().starts_with("EZ_"));
}
