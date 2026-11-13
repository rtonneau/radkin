use radkin::config::{Config, Section};

const F: &str = "tests/fixtures/config";

fn k_value(cfg: &Config, section: Section, reaction: &str) -> Option<toml::Value> {
    cfg.reactions[&section]
        .iter()
        .find(|r| r["reaction"].as_str() == Some(reaction))
        .and_then(|r| r.get("k_value").cloned())
}

#[test]
fn labarbe_config_loads_includes_and_beam() {
    let cfg = Config::load("data/config_Labarbe.toml").unwrap();
    let golden: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string("tests/golden/derive.json").unwrap())
            .unwrap();
    let n = &golden["labarbe"]["n_reactions"];
    assert_eq!(
        cfg.reactions[&Section::AcidBase].len() as u64,
        n["acid_base"].as_u64().unwrap()
    );
    assert_eq!(
        cfg.reactions[&Section::Equilibrium].len() as u64,
        n["equilibrium"].as_u64().unwrap()
    );
    let n_k = cfg.reactions[&Section::KReaction].len()
        + cfg.reactions[&Section::Arrhenius].len()
        + cfg.reactions[&Section::Constant].len();
    assert_eq!(n_k as u64, n["k"].as_u64().unwrap());
    let beams = radkin::physics::BeamCollection::from_config(&cfg.beam).unwrap();
    assert_eq!(beams.len(), 1);
    let (beam, let_) = beams.get("default").unwrap();
    assert_eq!(beam.dose_rate_peak, 1e7);
    assert_eq!(*let_, Some(0.201));
}

#[test]
fn override_replaces_reactions_in_place() {
    let cfg = Config::load(format!("{F}/override.toml")).unwrap();
    let k = &cfg.reactions[&Section::KReaction];
    assert_eq!(k.len(), 3);
    assert_eq!(
        k[0]["reaction"].as_str(),
        Some("B + A -> C"),
        "kept in place, string replaced"
    );
    assert_eq!(
        k_value(&cfg, Section::KReaction, "B + A -> C")
            .unwrap()
            .as_float(),
        Some(2e9)
    );
    assert_eq!(k[2]["reaction"].as_str(), Some("D -> E"));
    let eq = &cfg.reactions[&Section::Equilibrium][0];
    assert_eq!(eq["k_forward"].as_float(), Some(4.0));
    assert!(eq.get("k_backward").is_none());
    let ab = &cfg.reactions[&Section::AcidBase][0];
    assert_eq!(ab["k_forward"].as_float(), Some(0.2));
    assert_eq!(ab["pKa"].as_float(), Some(11.65));
    assert_eq!(cfg.temperature(), 37.0);
    assert_eq!(cfg.ph(), 7.0);
}

#[test]
fn nested_include_is_relative_to_including_file() {
    let cfg = Config::load(format!("{F}/sub/nested.toml")).unwrap();
    assert_eq!(cfg.initial["O2"].as_float(), Some(10e-6));
    assert!(
        cfg.fixed.contains_key("catalase"),
        "radiopyo/ prefix maps to data/"
    );
    assert_eq!(
        cfg.reactions[&Section::KReaction][0]["reaction"].as_str(),
        Some("A + B -> C")
    );
}

#[test]
fn include_cycle_is_an_error() {
    let err = Config::load(format!("{F}/cycle_a.toml"))
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("cycle_a.toml") && err.contains("cycle_b.toml"),
        "{err}"
    );
}

#[test]
fn missing_include_names_the_path() {
    let err = Config::load(format!("{F}/missing.toml"))
        .unwrap_err()
        .to_string();
    assert!(err.contains("does_not_exist.toml"), "{err}");
}

#[test]
fn from_table_matches_from_toml_str() {
    let text = std::fs::read_to_string(format!("{F}/base.toml")).unwrap();
    let base = std::path::Path::new(F);
    let table: toml::Table = toml::from_str(&text).unwrap();
    assert_eq!(
        Config::from_table(table, Some(base)).unwrap(),
        Config::from_toml_str(&text, Some(base)).unwrap()
    );
}
