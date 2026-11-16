//! Reactions of one particle in the IRT stage: first order, with a fixed
//! species, and with a background species.

use std::path::{Path, PathBuf};

use radkin::irt::{IrtOutput, IrtSimulation, run_irt_with};

/// Allowed distance in standard errors.
const K: f64 = 4.0;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_irt_bg_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

/// `n` molecules of `species` on a line, 1 um apart: no pair meets.
fn far_apart(dir: &Path, species: &str, n: usize) -> PathBuf {
    let mut csv = String::from("species,x,y,z\n");
    for i in 0..n {
        csv.push_str(&format!("{species},{},0,0\n", 1000.0 * i as f64));
    }
    let path = dir.join(format!("{species}.csv"));
    std::fs::write(&path, csv).unwrap();
    path
}

/// A run of `chemistry` (TOML text) with `[irt]` keys `irt` on `molecules`.
fn simulation(
    dir: &Path,
    name: &str,
    chemistry: &str,
    irt: &str,
    molecules: &Path,
) -> IrtSimulation {
    let text = format!(
        "{chemistry}\n[irt]\n{irt}\n[[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
        toml_path(molecules)
    );
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    IrtSimulation::from_toml(&path).unwrap()
}

fn index(out: &IrtOutput, label: &str) -> usize {
    out.species.iter().position(|s| s == label).unwrap()
}

/// `N0 exp(-rate t)` at every record time, for one trajectory of `n0`
/// molecules.
fn assert_decay(out: &IrtOutput, label: &str, n0: f64, rate: f64) {
    let s = index(out, label);
    for (t, row) in out.record_times.iter().zip(&out.totals_mean) {
        let p = (-rate * t).exp();
        let sd = (n0 * p * (1.0 - p)).sqrt().max(1e-9);
        assert!(
            (row[s] - n0 * p).abs() <= K * sd,
            "{label} at {t:e} s: {} against {} (sd {sd})",
            row[s],
            n0 * p
        );
    }
}

#[test]
fn irt_background_first_order_decay_is_exponential() {
    let dir = scratch("decay");
    let chemistry = "[[reactions.k_reaction]]\nreaction = \"HO2_minus -> H2O2\"\nk_value = 2e5";
    let n0 = 4_000;
    let mol = far_apart(&dir, "HO2_minus", n0);
    let irt = "t_end = 2e-5\nrecord_times = [1e-7, 1e-6, 3e-6, 1e-5]";
    let sim = simulation(&dir, "decay.toml", chemistry, irt, &mol);
    assert_eq!(sim.model.pairs.len(), 0);
    assert_eq!(sim.model.molecule.len(), 1);
    assert!((sim.model.molecule[0].rate - 2e5).abs() < 1e-9 * 2e5);
    let out = run_irt_with(&sim, 1, 1).unwrap();
    assert_decay(&out, "HO2_minus", n0 as f64, 2e5);
    // every molecule that reacted left one H2O2
    let (a, b) = (index(&out, "HO2_minus"), index(&out, "H2O2"));
    for row in &out.totals_mean {
        assert_eq!(row[a] + row[b], n0 as f64);
    }
    assert_eq!(
        out.stats.molecule_reactions as f64,
        out.totals_mean.last().unwrap()[b]
    );
    assert_eq!(out.stats.pair_reactions + out.stats.contact_reactions, 0);
    assert_eq!(out.stats.background_consumed.iter().sum::<u64>(), 0);
}

#[test]
fn irt_background_scavenging_by_a_background_species() {
    let dir = scratch("scavenging");
    // aerated water: 250 uM of O2, e_aq + O2 at 1.9e10 M-1 s-1
    let reaction =
        "[[reactions.k_reaction]]\nreaction = \"e_aq + O2 -> O2_minus\"\nk_value = 1.9e10";
    let rate = 1.9e10 * 2.5e-4;
    let n0 = 3_000;
    let mol = far_apart(&dir, "e_aq", n0);
    let irt = "t_end = 1e-6\nrecord_times = [1e-8, 1e-7, 3e-7]";

    // O2 as a species with an initial concentration: a background species
    let chemistry = format!("[concentrations.initial]\nO2 = 2.5e-4\n{reaction}");
    let sim = simulation(&dir, "initial.toml", &chemistry, irt, &mol);
    let o2 = sim.model.species.iter().position(|s| s == "O2").unwrap();
    assert_eq!(sim.model.background[o2], 2.5e-4);
    // one pair reaction (with O2 particles, none here) and its background twin
    assert_eq!(sim.model.pairs.len(), 1);
    assert_eq!(sim.model.molecule.len(), 1);
    assert_eq!(sim.model.molecule[0].background, Some(o2));
    assert!((sim.model.molecule[0].rate - rate).abs() < 1e-9 * rate);
    let out = run_irt_with(&sim, 1, 2).unwrap();
    assert_decay(&out, "e_aq", n0 as f64, rate);
    let last = out.totals_mean.last().unwrap();
    let reacted = n0 as f64 - last[index(&out, "e_aq")];
    assert_eq!(last[index(&out, "O2_minus")], reacted);
    // the O2 consumed is counted, and no O2 particle appears or disappears
    assert_eq!(out.stats.background_consumed[o2] as f64, reacted);
    assert_eq!(last[o2], 0.0);
    assert_eq!(out.stats.molecule_reactions as f64, reacted);

    // O2 as a fixed species: the same rate, nothing to count
    let chemistry = format!("[concentrations.fixed]\nO2 = 2.5e-4\n{reaction}");
    let sim = simulation(&dir, "fixed.toml", &chemistry, irt, &mol);
    assert_eq!(sim.model.pairs.len(), 0);
    assert_eq!(sim.model.molecule[0].background, None);
    assert!((sim.model.molecule[0].rate - rate).abs() < 1e-9 * rate);
    let out = run_irt_with(&sim, 1, 3).unwrap();
    assert_decay(&out, "e_aq", n0 as f64, rate);
    assert_eq!(out.stats.background_consumed.iter().sum::<u64>(), 0);
}

#[test]
fn irt_background_self_reaction_with_its_own_background() {
    let dir = scratch("self");
    // 2 OH_r -> H2O2 with k_obs 5.5e9 (the TOML holds 2 k_obs) and 10 uM of
    // OH_r in the background: a particle meets it at the rate 2 k_obs c
    let chemistry = "[concentrations.initial]\nOH_r = 1e-5\n\
                     [[reactions.k_reaction]]\nreaction = \"2 OH_r -> H2O2\"\nk_value = 1.1e10";
    let rate = 2.0 * 5.5e9 * 1e-5;
    let n0 = 2_000;
    let mol = far_apart(&dir, "OH_r", n0);
    // far apart, so the pair reaction between particles plays no role
    let irt = "t_end = 2e-5\ncutoff_probability = 0.5\nrecord_times = [1e-6, 5e-6, 1e-5]";
    let sim = simulation(&dir, "self.toml", chemistry, irt, &mol);
    assert_eq!(sim.model.molecule.len(), 1);
    assert!((sim.model.molecule[0].rate - rate).abs() < 1e-9 * rate);
    let out = run_irt_with(&sim, 1, 4).unwrap();
    assert_eq!(out.stats.pair_reactions, 0);
    assert_decay(&out, "OH_r", n0 as f64, rate);
    let oh = index(&out, "OH_r");
    assert_eq!(
        out.stats.background_consumed[oh] as f64,
        n0 as f64 - out.totals_mean.last().unwrap()[oh]
    );
}

#[test]
fn irt_background_switch_turns_every_single_molecule_reaction_off() {
    let dir = scratch("switch");
    let chemistry = "[concentrations.initial]\nO2 = 2.5e-4\n\
                     [[reactions.k_reaction]]\nreaction = \"e_aq + O2 -> O2_minus\"\nk_value = 1.9e10\n\
                     [[reactions.k_reaction]]\nreaction = \"e_aq -> H_r\"\nk_value = 1e7";
    let mol = far_apart(&dir, "e_aq", 500);
    let sim = simulation(
        &dir,
        "off.toml",
        chemistry,
        "t_end = 1e-5\nbackground_reactions = false",
        &mol,
    );
    let out = run_irt_with(&sim, 2, 5).unwrap();
    assert_eq!(out.stats.molecule_reactions, 0);
    let e = index(&out, "e_aq");
    assert!(out.totals_mean.iter().all(|row| row[e] == 500.0));
    // on (the default): every e_aq is gone by 10 us, by one way or the other
    let sim = simulation(&dir, "on.toml", chemistry, "t_end = 1e-5", &mol);
    let out = run_irt_with(&sim, 2, 5).unwrap();
    let last = out.totals_mean.last().unwrap();
    assert_eq!(last[e], 0.0);
    assert_eq!(out.stats.molecule_reactions, 1000);
    // the two ways share the molecules as their rates: 4.75e6 against 1e7
    let by_o2 = last[index(&out, "O2_minus")] / 500.0;
    let p = 4.75e6 / (4.75e6 + 1e7);
    assert!(
        (by_o2 - p).abs() < K * (p * (1.0 - p) / 1000.0).sqrt(),
        "{by_o2} against {p}"
    );
}

#[test]
fn irt_background_ph7_water_neutralises_the_track_ions() {
    // Ten Geant4 events in free space, with and without the background
    // reactions of data/chem_g4.toml (water and the pH 7 ions).
    let dir = scratch("ph7");
    let molecules = root().join("tests/data/irt/g4_events_molecules.csv");
    let run = |name: &str, on: bool| {
        let text = format!(
            "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nt_end = 1e-5\nstart_contact = false\n\
             background_reactions = {on}\nrecord_times = [1e-7]\n\
             [[irt.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [1, 4, 5]\n",
            toml_path(&molecules)
        );
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        let sim = IrtSimulation::from_toml(&path).unwrap();
        run_irt_with(&sim, 48, 6).unwrap()
    };
    let off = run("off.toml", false);
    let on = run("on.toml", true);
    assert_eq!(off.stats.molecule_reactions, 0);
    assert!(on.stats.molecule_reactions > 0);
    let h3o = index(&on, "H3O_plus");
    let oh_minus = index(&on, "OH_minus");
    // at 100 ns the background has not acted yet (rates up to 1.1e4 s-1)
    let sem = |o: &IrtOutput, k: usize, s: usize| o.totals_std[k][s] / 48f64.sqrt();
    let diff = on.totals_mean[1][h3o] - off.totals_mean[1][h3o];
    assert!(
        diff.abs() <= K * (sem(&on, 1, h3o).powi(2) + sem(&off, 1, h3o).powi(2)).sqrt() + 0.5,
        "{diff}"
    );
    // at 10 us about 1 - exp(-1.13e11 x 1e-7 x 1e-5) = 11 % of the H3O+ left
    // at 100 ns has met a background OH-
    let ratio = on.totals_mean[2][h3o] / off.totals_mean[2][h3o];
    assert!((0.84..0.95).contains(&ratio), "{ratio}");
    assert!(on.stats.background_consumed[oh_minus] > 0);
    assert!(on.stats.background_consumed[h3o] > 0);
    // the background is not made of particles
    let meta = radkin::irt::output::irt_meta_json(
        &IrtSimulation::from_toml(dir.join("on.toml")).unwrap(),
        &on,
    )
    .unwrap();
    let meta: serde_json::Value = serde_json::from_str(&meta).unwrap();
    assert!((meta["background"]["H3O_plus"].as_f64().unwrap() - 1e-7).abs() < 1e-20);
    assert!(meta["background_consumed"]["OH_minus"].as_u64().unwrap() > 0);
    assert_eq!(meta["background_reactions"], true);
}
