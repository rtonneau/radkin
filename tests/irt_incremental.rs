//! Incremental resampling of the IRT stage against full resampling, on the
//! same seeds.

use std::path::{Path, PathBuf};

use radkin::irt::{IrtOutput, IrtSimulation, run_irt_with};
use radkin::ssa::SsaRng;

/// Allowed distance in standard errors.
const K: f64 = 4.0;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_irt_incr_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

/// The same run with the two rules: `text` holds `RULE` where the value of
/// `resampling` goes.
fn both(dir: &Path, text: &str, trajectories: usize, seed: u64) -> (IrtOutput, IrtOutput) {
    let run = |rule: &str, seed: u64| {
        let path = dir.join(format!("{rule}.toml"));
        std::fs::write(&path, text.replace("RULE", rule)).unwrap();
        let sim = IrtSimulation::from_toml(&path).unwrap();
        run_irt_with(&sim, trajectories, seed).unwrap()
    };
    (run("full", seed), run("incremental", seed + 1))
}

/// Every species at every record within `K` standard errors of the
/// difference (plus `floor` molecules); returns the largest |z|.
fn assert_same(full: &IrtOutput, incremental: &IrtOutput, trajectories: usize, floor: f64) -> f64 {
    assert_eq!(full.record_times, incremental.record_times);
    let n = trajectories as f64;
    let mut worst: f64 = 0.0;
    for (k, t) in full.record_times.iter().enumerate() {
        for (s, label) in full.species.iter().enumerate() {
            let (a, b) = (full.totals_mean[k][s], incremental.totals_mean[k][s]);
            let sd =
                ((full.totals_std[k][s].powi(2) + incremental.totals_std[k][s].powi(2)) / n).sqrt();
            assert!(
                (a - b).abs() <= K * sd + floor,
                "{label} at {t:e} s: full {a}, incremental {b} (sd {sd})"
            );
            if sd > 0.0 {
                worst = worst.max((a - b).abs() / sd);
            }
        }
    }
    worst
}

fn molecules_file() -> String {
    toml_path(&root().join("tests/data/irt/g4_events_molecules.csv"))
}

const HEAD: &str = "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nresampling = \"RULE\"\nstart_contact = false\n";

#[test]
fn irt_incremental_agrees_with_full_on_single_tracks() {
    let dir = scratch("single");
    let trajectories = 48;
    let mut worst: f64 = 0.0;
    for event in [2, 6, 9] {
        let text = format!(
            "{HEAD}background_reactions = false\nt_end = 1e-6\nrecord_times = [1e-11, 1e-10, 1e-9, 1e-8, 1e-7]\n\
             [[irt.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [{event}]\n",
            molecules_file()
        );
        let (full, incremental) = both(&dir, &text, trajectories, 70 + event);
        assert!(incremental.stats.pair_reactions > 0);
        worst = worst.max(assert_same(&full, &incremental, trajectories, 0.3));
    }
    println!("single tracks: largest |z| {worst:.2}");
}

#[test]
fn irt_incremental_agrees_with_full_on_overlapping_tracks() {
    let dir = scratch("overlap");
    // five tracks through the same 300 nm wide region at the same time, with
    // the background reactions of water at pH 7
    let mut seeds = String::new();
    for (k, event) in [0, 1, 4, 5, 8].iter().enumerate() {
        seeds.push_str(&format!(
            "[[irt.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [{event}]\noffset_nm = [{}, {}, 0]\n",
            molecules_file(),
            60.0 * k as f64 - 120.0,
            -45.0 * k as f64 + 90.0
        ));
    }
    let text = format!("{HEAD}t_end = 1e-6\nrecord_times = [1e-10, 1e-9, 1e-8, 1e-7]\n{seeds}");
    let trajectories = 24;
    let (full, incremental) = both(&dir, &text, trajectories, 80);
    let worst = assert_same(&full, &incremental, trajectories, 0.5);
    println!("overlapping tracks: largest |z| {worst:.2}");
    // the tracks do overlap: more reactions than the five alone would have
    assert!(incremental.stats.molecule_reactions > 0);
    assert!(
        incremental.stats.pair_reactions > 24 * 60,
        "{:?}",
        incremental.stats
    );
}

#[test]
fn irt_incremental_agrees_with_full_between_reflective_walls() {
    let dir = scratch("walls");
    // A + B -> C, partially diffusion-controlled, in a 40 nm box
    let mut rng = SsaRng::new(90, 0);
    let mut csv = String::from("species,x,y,z\n");
    for i in 0..40 {
        let name = if i < 20 { "O2" } else { "O3" };
        let p: [f64; 3] = std::array::from_fn(|_| 40.0 * rng.uniform());
        csv.push_str(&format!("{name},{},{},{}\n", p[0], p[1], p[2]));
    }
    let mol = dir.join("box.csv");
    std::fs::write(&mol, csv).unwrap();
    let text = format!(
        "[[reactions.k_reaction]]\nreaction = \"O2 + O3 -> O_atom\"\nk_value = 3e9\nirt_type = \"partial\"\n\
         [irt]\nresampling = \"RULE\"\nstart_contact = false\nboundary = \"reflective\"\nt_end = 1.6e-6\n\
         record_times = [2e-7, 5e-7, 1e-6, 1.5e-6]\n[irt.box]\nlo_nm = 0\nhi_nm = 40\n\
         [irt.species.O2]\nradius_nm = 0.25\ndiffusion = 2e-9\n[irt.species.O3]\nradius_nm = 0.25\ndiffusion = 2e-9\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
        toml_path(&mol)
    );
    let trajectories = 150;
    let (full, incremental) = both(&dir, &text, trajectories, 91);
    let worst = assert_same(&full, &incremental, trajectories, 0.1);
    println!("reflective box: largest |z| {worst:.2}");
    // both rules draw every pair again at the look-ahead limit of the box
    assert!(full.stats.box_moves > 0 && incremental.stats.box_moves > 0);
    assert!(incremental.stats.image_reactions > 0);
    // about half of the molecules have reacted
    let a = full.species.iter().position(|s| s == "O2").unwrap();
    let left = incremental.totals_mean.last().unwrap()[a];
    assert!((4.0..14.0).contains(&left), "{left}");
}

#[test]
fn irt_incremental_single_molecule_reactions_and_contact_rule() {
    let dir = scratch("alone");
    // first-order decay of molecules far apart
    let mut csv = String::from("species,x,y,z\n");
    for i in 0..3000 {
        csv.push_str(&format!("HO2_minus,{},0,0\n", 1000.0 * i as f64));
    }
    let mol = dir.join("decay.csv");
    std::fs::write(&mol, csv).unwrap();
    let text = format!(
        "[[reactions.k_reaction]]\nreaction = \"HO2_minus -> H2O2\"\nk_value = 2e5\n\
         [irt]\nresampling = \"incremental\"\nt_end = 2e-5\nrecord_times = [1e-6, 3e-6, 1e-5]\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
        toml_path(&mol)
    );
    let path = dir.join("decay.toml");
    std::fs::write(&path, text).unwrap();
    let out = run_irt_with(&IrtSimulation::from_toml(&path).unwrap(), 1, 95).unwrap();
    let s = out.species.iter().position(|x| x == "HO2_minus").unwrap();
    for (t, row) in out.record_times.iter().zip(&out.totals_mean) {
        let p = (-2e5 * t).exp();
        let sd = (3000.0 * p * (1.0 - p)).sqrt().max(1e-9);
        assert!((row[s] - 3000.0 * p).abs() <= K * sd, "{t:e}: {}", row[s]);
    }
    assert_eq!(
        out.stats.molecule_reactions as f64,
        3000.0 - out.totals_mean.last().unwrap()[s]
    );
    // contact rule at the injection: H3O+ and OH- 0.1 nm apart
    let mol = dir.join("pair.csv");
    std::fs::write(&mol, "species,x,y,z\nH3O_plus,0,0,0\nOH_minus,0.1,0,0\n").unwrap();
    let text = format!(
        "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nresampling = \"incremental\"\nt_end = 1e-9\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 1e-12\n",
        toml_path(&mol)
    );
    let path = dir.join("pair.toml");
    std::fs::write(&path, &text).unwrap();
    let out = run_irt_with(&IrtSimulation::from_toml(&path).unwrap(), 8, 96).unwrap();
    let h = out.species.iter().position(|x| x == "H3O_plus").unwrap();
    assert_eq!(out.totals_mean[0][h], 1.0);
    assert_eq!(out.totals_mean[1][h], 0.0);
    assert_eq!(out.stats.contact_reactions, 8);
    // without the rule the pair is looked at again after it has moved
    let path = dir.join("pair_off.toml");
    std::fs::write(
        &path,
        text.replace("t_end = 1e-9", "t_end = 1e-9\nstart_contact = false"),
    )
    .unwrap();
    let out = run_irt_with(&IrtSimulation::from_toml(&path).unwrap(), 200, 97).unwrap();
    assert_eq!(out.stats.contact_reactions, 0);
    assert!(out.stats.moves >= 200, "{:?}", out.stats);
    let left = out.totals_mean.last().unwrap()[h];
    assert!((0.3..0.8).contains(&left), "{left}");
}

#[test]
fn irt_incremental_same_seed_same_result_and_meta() {
    let dir = scratch("repro");
    let text = format!(
        "{HEAD}t_end = 1e-7\n[[irt.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [7]\n",
        molecules_file()
    )
    .replace("RULE", "incremental");
    let path = dir.join("e.toml");
    std::fs::write(&path, text).unwrap();
    let sim = IrtSimulation::from_toml(&path).unwrap();
    let a = run_irt_with(&sim, 5, 3).unwrap();
    let b = run_irt_with(&sim, 5, 3).unwrap();
    assert_eq!(a.totals_mean, b.totals_mean);
    assert_eq!(a.stats, b.stats);
    assert_ne!(a.totals_mean, run_irt_with(&sim, 5, 4).unwrap().totals_mean);
    let meta = radkin::irt::output::irt_meta_json(&sim, &a).unwrap();
    let meta: serde_json::Value = serde_json::from_str(&meta).unwrap();
    assert_eq!(meta["resampling"], "incremental");
    // an unknown rule is refused
    let bad = "[irt]\nt_end = 1e-6\nresampling = \"lazy\"";
    let e = radkin::irt::IrtConfig::from_toml_str(bad)
        .unwrap_err()
        .to_string();
    assert!(e.contains("\"full\" or \"incremental\""), "{e}");
}
