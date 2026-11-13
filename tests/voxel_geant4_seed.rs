//! A real Geant4-DNA 11.4.3 option3 event (`tests/data/geant4_option3_event.csv`) seeded into
//! the chemistry of `data/voxel_example.toml`.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::consts::data_dir;
use radkin::voxel::{VoxelConfig, load_seeds, seed::parse_molecules};

const FIXTURE_ROWS: usize = 9672;
/// Geant4 names present in the fixture.
const GEANT4_SPECIES: [&str; 12] = [
    "e_aq^-1", "°OH^0", "H3O^1", "H^0", "H_2^0", "H2O2^0", "OH^-1", "HO_2°^0", "HO_2^-1", "O_2^0",
    "O^-1", "°O^0",
];

fn voxel_example_env() -> SimEnv {
    let cfg = Config::load(data_dir().join("voxel_example.toml")).unwrap();
    SimEnv::from_config(&cfg).unwrap()
}

fn fixture_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("geant4_option3_event.csv")
}

/// Run TOML: a 400 nm cube, the fixture shifted so its first `°O^0` sits at the cube centre.
/// Everything farther away (the track is ~70 um long) is clipped on purpose.
fn run_toml(clip: bool) -> String {
    format!(
        r#"
        [voxel]
        clip = {clip}
        [voxel.grid]
        dims = [40, 40, 40]
        spacing_nm = 10.0
        [voxel.geant4]
        ignore = ["H2O^0"]
        [[voxel.seed]]
        file = {file:?}
        time = 1.001e-12
        offset_nm = [73295.37512, 222933.8576, 243231.8982]
        "#,
        file = fixture_path().to_str().unwrap()
    )
}

#[test]
fn fixture_holds_every_geant4_species() {
    let mols = parse_molecules(&std::fs::read_to_string(fixture_path()).unwrap()).unwrap();
    assert_eq!(mols.len(), FIXTURE_ROWS);
    for name in GEANT4_SPECIES {
        assert!(mols.iter().any(|m| m.species == name), "{name} missing");
    }
    let mut names: Vec<&str> = mols.iter().map(|m| m.species.as_str()).collect();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), GEANT4_SPECIES.len());
}

#[test]
fn geant4_event_seeds_into_voxel_example_chemistry() {
    let env = voxel_example_env();
    let dir = std::env::temp_dir().join(format!("radkin-voxel-g4-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let toml = dir.join("run.toml");
    std::fs::write(&toml, run_toml(true)).unwrap();
    let vc = VoxelConfig::from_toml_file(&toml).unwrap();
    // Every row either resolves to a dynamic species (an unknown name would be an error, also
    // for the clipped rows) and is placed, or is clipped on purpose.
    let seeds = load_seeds(&vc.seeds, &vc.grid, &env.species, &vc.policy).unwrap();
    let rep = &seeds[0].report;
    assert_eq!(rep.placed + rep.clipped + rep.ignored, FIXTURE_ROWS);
    assert_eq!(rep.ignored, 0);
    assert!(rep.placed >= 300, "placed {}", rep.placed);
    assert!(rep.clipped > 0);
    assert_eq!(rep.per_species.values().sum::<usize>(), rep.placed);
    assert!(rep.per_species["O_atom"] >= 1);
    // The same seed without clipping is refused.
    std::fs::write(&toml, run_toml(false)).unwrap();
    let strict = VoxelConfig::from_toml_file(&toml).unwrap();
    let err = load_seeds(&strict.seeds, &strict.grid, &env.species, &strict.policy).unwrap_err();
    assert!(err.to_string().contains("outside the grid"), "{err}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn voxel_example_has_the_option3_o_atom_reactions() {
    let env = voxel_example_env();
    assert!(env.species.is_dyn("O_atom"));
    let strings = env.reaction_strings();
    let k = env.k_values().unwrap();
    let with_o: Vec<(&String, f64)> = strings
        .iter()
        .zip(&k)
        .filter(|(s, _)| s.contains("O_atom"))
        .map(|(s, k)| (s, *k))
        .collect();
    assert_eq!(with_o.len(), 10, "{with_o:?}");
    // O_atom + H2O -> 2 OH_r: pseudo-first-order 1.9e3 /s once [H2O] is multiplied in.
    let (_, k_h2o) = with_o
        .iter()
        .find(|(s, _)| s.contains("H2O") && !s.contains("H2O2"))
        .unwrap();
    let h2o = env.species.get("H2O").unwrap().initial_cc(env.temperature);
    assert!((k_h2o * h2o / 1.9e3 - 1.0).abs() < 1e-3, "{}", k_h2o * h2o);
}
