//! Reaction model of the IRT stage against the reaction table Geant4-DNA
//! 11.4.3 (option3) dumped at run time.

use std::collections::HashMap;
use std::path::PathBuf;

use radkin::chemistry::env::SimEnv;
use radkin::chemistry::reactions::IrtType;
use radkin::config::Config;
use radkin::irt::{IrtModel, SpeciesProps, SpeciesTable, TypeSource};
use radkin::matrices::SimMatrices;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Geant4 user IDs -> radkin labels.
const NAMES: [(&str, &str); 15] = [
    ("e_aq", "e_aq"),
    ("°OH", "OH_r"),
    ("H3Op", "H3O_plus"),
    ("H", "H_r"),
    ("H2", "H2"),
    ("H2O2", "H2O2"),
    ("OHm", "OH_minus"),
    ("HO2°", "HO2"),
    ("HO2m", "HO2_minus"),
    ("O2", "O2"),
    ("O2m", "O2_minus"),
    ("O3", "O3"),
    ("O3m", "O3_minus"),
    ("Om", "O_minus"),
    ("Oxy", "O_atom"),
];

struct DumpSpecies {
    diffusion: f64,
    charge: i32,
    radius_nm: f64,
}

struct DumpReaction {
    id: usize,
    a: String,
    b: String,
    kind: u8,
    k_obs: f64,
    r_reaction_nm: f64,
    r_eff_nm: f64,
    products: Vec<String>,
}

fn dump() -> (HashMap<String, DumpSpecies>, Vec<DumpReaction>) {
    let text =
        std::fs::read_to_string(root().join("tests/data/irt/g4_chemistry_dump.txt")).unwrap();
    let names: HashMap<&str, &str> = NAMES.into_iter().collect();
    let mut species = HashMap::new();
    let mut reactions = Vec::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        match p[0] {
            "species" => {
                if let Some(label) = names.get(p[2]) {
                    species.insert(
                        label.to_string(),
                        DumpSpecies {
                            diffusion: p[3].parse().unwrap(),
                            charge: p[4].parse().unwrap(),
                            radius_nm: p[5].parse().unwrap(),
                        },
                    );
                }
            }
            "reaction" => {
                // the others are scavenger reactions: NoneM, H2O(B), H3Op(B), OHm(B)
                let (Some(a), Some(b)) = (names.get(p[2]), names.get(p[3])) else {
                    continue;
                };
                reactions.push(DumpReaction {
                    id: p[1].parse().unwrap(),
                    a: a.to_string(),
                    b: b.to_string(),
                    kind: p[4].parse().unwrap(),
                    k_obs: p[6].parse().unwrap(),
                    r_reaction_nm: p[8].parse().unwrap(),
                    r_eff_nm: p[9].parse().unwrap(),
                    products: if p[10] == "-" {
                        Vec::new()
                    } else {
                        p[10].split(';').map(|x| names[x].to_string()).collect()
                    },
                });
            }
            _ => {}
        }
    }
    (species, reactions)
}

fn g4_text() -> String {
    std::fs::read_to_string(root().join("data/chem_g4.toml")).unwrap()
}

fn model_of(toml: &str, table: &SpeciesTable) -> IrtModel {
    let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    IrtModel::new(&env, &mat, table).unwrap()
}

fn close(a: f64, b: f64, rtol: f64) -> bool {
    (a - b).abs() <= rtol * b.abs()
}

#[test]
fn irt_model_species_table_equals_the_geant4_dump() {
    let (species, _) = dump();
    let table = SpeciesTable::water();
    assert_eq!(table.len(), species.len());
    for (label, d) in &species {
        let p = table
            .get(label)
            .unwrap_or_else(|| panic!("{label} missing"));
        assert_eq!(p.diffusion, d.diffusion, "{label}");
        assert_eq!(p.charge, d.charge, "{label}");
        assert!(close(p.radius, d.radius_nm * 1e-9, 1e-15), "{label}");
    }
}

#[test]
fn irt_model_geant4_network_has_geant4_types_and_radii() {
    let (_, reactions) = dump();
    assert_eq!(reactions.len(), 52);
    let model = model_of(&g4_text(), &SpeciesTable::water());
    assert_eq!(model.pairs.len(), 52);
    assert!(model.warnings.is_empty(), "{:?}", model.warnings);
    let idx = |l: &str| model.species.iter().position(|s| s == l).unwrap();
    for r in &reactions {
        let ids = model.pairs_of(idx(&r.a), idx(&r.b));
        assert_eq!(ids.len(), 1, "G4 reaction {}", r.id);
        let p = &model.pairs[ids[0] as usize];
        let what = format!("G4 reaction {} ({} + {})", r.id, r.a, r.b);
        assert_eq!(p.source, TypeSource::Network, "{what}");
        let kind = if r.kind == 0 {
            IrtType::Full
        } else {
            IrtType::Partial
        };
        assert_eq!(p.kind, kind, "{what}");
        assert!(close(p.k_obs, r.k_obs, 1e-12), "{what}: k_obs {}", p.k_obs);
        assert!(
            close(p.r_eff, r.r_eff_nm * 1e-9, 1e-9),
            "{what}: r_eff {} nm against {}",
            p.r_eff * 1e9,
            r.r_eff_nm
        );
        assert!(
            close(p.r_reaction, r.r_reaction_nm * 1e-9, 1e-9),
            "{what}: r_reaction {} nm against {}",
            p.r_reaction * 1e9,
            r.r_reaction_nm
        );
        let mut got: Vec<&str> = p
            .products
            .iter()
            .map(|&s| model.species[s].as_str())
            .collect();
        let mut want: Vec<&str> = r.products.iter().map(|s| s.as_str()).collect();
        got.sort_unstable();
        want.sort_unstable();
        assert_eq!(got, want, "{what}");
        if kind == IrtType::Partial {
            assert!(p.k_obs < p.k_dif, "{what}");
            assert!(p.k_act > 0.0 && p.k_act.is_finite(), "{what}");
            assert!(p.p_contact > 0.0 && p.p_contact < 1.0, "{what}");
            // k_obs is the series of the diffusion and activation steps
            assert!(
                close(p.k_dif * p.k_act / (p.k_dif + p.k_act), p.k_obs, 1e-12),
                "{what}"
            );
        } else {
            assert_eq!(p.p_contact, 1.0, "{what}");
            assert_eq!(p.r_onsager, 0.0, "{what}");
        }
    }
}

#[test]
fn irt_model_self_reactions_use_the_event_rate_constant() {
    // "2 OH_r -> H2O2" holds 2 k_obs in the TOML; Geant4's k_obs is 5.5e9.
    let model = model_of(&g4_text(), &SpeciesTable::water());
    let oh = model.species.iter().position(|s| s == "OH_r").unwrap();
    let p = &model.pairs[model.pairs_of(oh, oh)[0] as usize];
    assert!(close(p.k_obs, 5.5e9, 1e-12));
    // identical reactants: half the diffusion rate of two distinct species
    let d = model.diffusion[oh];
    let k_dif = 4.0 * std::f64::consts::PI * (2.0 * d) * p.r_eff * 6.02214076e23 * 1e3 / 2.0;
    assert!(close(p.k_dif, k_dif, 1e-12));
}

/// The rule `k_obs >= k_dif` on the Geant4 network, against Geant4's types.
#[test]
fn irt_model_rate_constant_rule_against_geant4_types() {
    let untyped: String = g4_text()
        .lines()
        .filter(|l| !l.starts_with("irt_type"))
        .collect::<Vec<_>>()
        .join("\n");
    let typed = model_of(&g4_text(), &SpeciesTable::water());
    let ruled = model_of(&untyped, &SpeciesTable::water());
    assert_eq!(ruled.pairs.len(), 52);
    let mut same = 0;
    let mut full_in_g4 = Vec::new();
    let mut partial_in_g4 = Vec::new();
    for (t, r) in typed.pairs.iter().zip(&ruled.pairs) {
        assert_eq!(r.source, TypeSource::Rule);
        let name = format!("{} + {}", typed.species[t.a], typed.species[t.b]);
        match (t.kind, r.kind) {
            (a, b) if a == b => same += 1,
            (IrtType::Full, _) => full_in_g4.push(name),
            (IrtType::Partial, _) => partial_in_g4.push(name),
        }
    }
    println!("rule agrees with Geant4 on {same}/52 pair reactions");
    println!("full in Geant4, partial by the rule: {full_in_g4:?}");
    println!("partial in Geant4, full by the rule: {partial_in_g4:?}");
    // Geant4 never types "partial" a reaction at or above its diffusion limit.
    assert!(partial_in_g4.is_empty());
    assert_eq!(same, 47);
    assert_eq!(
        full_in_g4,
        [
            "H_r + H_r",
            "H_r + e_aq",
            "H_r + O_atom",
            "H_r + O_minus",
            "e_aq + e_aq"
        ]
    );
}

#[test]
fn irt_model_species_missing_from_the_table_falls_back_to_full() {
    let toml = r#"
[[reactions.k_reaction]]
reaction = "e_aq + Cl_r -> Cl_minus"
k_value = 1e10
[[reactions.k_reaction]]
reaction = "e_aq + OH_r -> OH_minus"
k_value = 2.95e10
"#;
    let mut table = SpeciesTable::water();
    let model = model_of(toml, &table);
    assert_eq!(model.pairs.len(), 2);
    let cl = &model.pairs[0];
    assert_eq!((cl.kind, cl.source), (IrtType::Full, TypeSource::Fallback));
    assert_eq!(model.pairs[1].source, TypeSource::Rule);
    assert_eq!(model.warnings.len(), 1);
    assert!(model.warnings[0].contains("Cl_r"), "{}", model.warnings[0]);

    // with an entry the rule applies
    table.set(
        "Cl_r",
        SpeciesProps {
            radius: 0.2e-9,
            charge: 0,
            diffusion: 2.0e-9,
        },
    );
    let model = model_of(toml, &table);
    assert_eq!(model.pairs[0].source, TypeSource::Rule);
    assert!(model.warnings.is_empty());
}

#[test]
fn irt_model_molecule_reactions_first_order_water_and_background() {
    let model = model_of(&g4_text(), &SpeciesTable::water());
    let idx = |l: &str| model.species.iter().position(|s| s == l).unwrap();
    let rate_of = |species: &str, background: Option<&str>, n_products: usize| -> Vec<f64> {
        model
            .molecule_of(idx(species))
            .iter()
            .map(|&i| &model.molecule[i as usize])
            .filter(|m| m.background == background.map(idx) && m.products.len() == n_products)
            .map(|m| m.rate)
            .collect()
    };
    // O3_minus -> H3O_plus + O2_minus: first order, 2660 s-1 in Geant4
    assert!(close(rate_of("O3_minus", None, 2)[0], 2660.0, 1e-9));
    // e_aq + H2O -> H_r + OH_minus: 15.8 s-1 at Geant4's [H2O] of 55.93 M; radkin's
    // water concentration at 25 degC is 55.9 M
    let r = rate_of("e_aq", None, 2)[0];
    assert!(close(r, 15.8, 2e-3), "{r}");
    // pH 7 background: e_aq + H3O_plus(background) at k 2.11e10 x 1e-7 M
    let r = rate_of("e_aq", Some("H3O_plus"), 1);
    assert_eq!(r.len(), 1);
    assert!(close(r[0], 2.11e10 * 1e-7, 1e-9), "{}", r[0]);
    // water autoionisation has no particle reactant
    assert_eq!(model.skipped.len(), 1);
    assert!(model.skipped[0].1.contains("H2O"), "{}", model.skipped[0].1);
}

#[test]
fn irt_model_unsupported_networks_are_refused() {
    let build = |toml: &str| {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        IrtModel::new(&env, &mat, &SpeciesTable::water())
            .unwrap_err()
            .to_string()
    };
    let e =
        build("[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r + H_r -> H2O2\"\nk_value = 1e10");
    assert!(e.contains("more than two particle reactants"), "{e}");
    let e = build("[[reactions.k_reaction]]\nreaction = \"A_x + B_x -> C_x\"\nk_value = 1e10");
    assert!(e.contains("do not diffuse"), "{e}");
    let bad = Config::from_toml_str(
        "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = 1e10\nirt_type = \"fast\"",
        None,
    )
    .and_then(|c| SimEnv::from_config(&c))
    .unwrap_err()
    .to_string();
    assert!(bad.contains("irt_type"), "{bad}");
}
