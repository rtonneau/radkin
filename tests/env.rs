use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::gvalues::GValuesBoscolo;
use serde_json::Value;

fn golden() -> Value {
    serde_json::from_str(&std::fs::read_to_string("tests/golden/derive.json").unwrap()).unwrap()
}

fn config(name: &str) -> &'static str {
    match name {
        "labarbe" => "data/config_Labarbe.toml",
        "pure_water" => "data/advanced_water_radiolysis.toml",
        _ => unreachable!(),
    }
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

fn close(a: f64, b: f64, rtol: f64) -> bool {
    a == b || (a - b).abs() <= rtol * b.abs()
}

#[test]
fn species_reactions_and_labels_match_radiopyo() {
    let g = golden();
    for name in ["labarbe", "pure_water"] {
        let env = SimEnv::from_config(&Config::load(config(name)).unwrap()).unwrap();
        let e = &g[name];
        assert_eq!(
            env.species.sorted_dyn_labels(),
            strings(&e["species_order"]),
            "{name} dyn order"
        );
        assert_eq!(
            env.species.sorted_cst_labels(),
            strings(&e["cst_species_order"]),
            "{name} cst order"
        );
        assert_eq!(
            env.reaction_strings(),
            strings(&e["reaction_str"]),
            "{name} str(reaction)"
        );
        assert_eq!(
            env.reaction_labels(),
            strings(&e["reaction_label"]),
            "{name} labels"
        );
        let k = env.k_values().unwrap();
        for (i, (a, b)) in k.iter().zip(e["k_values"].as_array().unwrap()).enumerate() {
            assert!(
                close(*a, b.as_f64().unwrap(), 1e-13),
                "{name} k[{i}] {a} vs {b}"
            );
        }
        let y0 = env.initial_values();
        for (a, b) in y0.iter().zip(e["initial_values"].as_array().unwrap()) {
            assert!(
                close(*a, b.as_f64().unwrap(), 1e-14),
                "{name} y0 {a} vs {b}"
            );
        }
    }
}

#[test]
fn let_sets_boscolo_g_values() {
    let g = golden();
    let gv = GValuesBoscolo::load(&radkin::consts::data_dir()).unwrap();
    for name in ["labarbe", "pure_water"] {
        let mut env = SimEnv::from_config(&Config::load(config(name)).unwrap()).unwrap();
        env.apply_let(&gv, 0.201, None).unwrap();
        for (label, expected) in g[name]["g_values"].as_object().unwrap() {
            let r = env
                .reactions
                .radiolytic
                .iter()
                .find(|r| &r.label() == label)
                .unwrap();
            assert!(
                close(r.ge(), expected.as_f64().unwrap(), 1e-12),
                "{name} {label}: {} vs {expected}",
                r.ge()
            );
        }
    }
}

#[test]
fn nk_convention_and_arrhenius() {
    let cfg = Config::from_toml_str(
        "[parameters]\n[[reactions.k_reaction]]\nreaction = \"2 OH_r -> H2O2\"\nk_value = 1.1e10\n\
         [[reactions.arrhenius]]\nreaction = \"H_r + OH_r -> H2O\"\nA = 1e10\nEa = 20\n",
        None,
    )
    .unwrap();
    let env = SimEnv::from_config(&cfg).unwrap();
    let k = env.k_values().unwrap();
    assert!(close(k[0], 5.5e9, 1e-14), "{}", k[0]);
    assert!(close(k[1], 3.13e6, 1e-2), "{}", k[1]);
    assert_eq!(env.reaction_strings()[0], "2OH_r -> H2O2");
}
