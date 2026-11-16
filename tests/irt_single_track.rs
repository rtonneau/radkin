//! IRT stage on single tracks: bookkeeping, the outflow box, the CLI, and
//! Geant4-DNA IRT_syn on ten reference events (`tests/data/irt/`).

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use radkin::irt::{IrtSimulation, SpeciesTable, run_irt_with};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_irt_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

/// A config of the Geant4 chemistry with the `[irt]` text `irt`.
fn config(dir: &Path, name: &str, irt: &str) -> PathBuf {
    let path = dir.join(name);
    let text = format!("[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n\n{irt}\n");
    std::fs::write(&path, text).unwrap();
    path
}

/// The ten Geant4 reference events.
struct Reference {
    molecules: PathBuf,
    /// Deposited energy [eV] by event.
    edep: HashMap<i64, f64>,
    /// `(event, species, time index)` -> molecules; times of `TIMES`.
    counts: HashMap<(i64, String, usize), f64>,
}

/// Times of the reference counts [s]; the first one is the start.
const TIMES: [f64; 7] = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6];
/// Time of the seed molecules [s].
const T_SEED: f64 = 1.001e-12;

fn reference() -> Reference {
    let dir = root().join("tests/data/irt");
    let read = |name: &str| std::fs::read_to_string(dir.join(name)).unwrap();
    let edep = read("g4_events_edep.csv")
        .lines()
        .skip(1)
        .map(|l| {
            let (e, v) = l.split_once(',').unwrap();
            (e.parse().unwrap(), v.parse().unwrap())
        })
        .collect();
    let mut counts = HashMap::new();
    for l in read("g4_events_gt.csv").lines().skip(1) {
        let f: Vec<&str> = l.split(',').collect();
        let t: f64 = f[2].parse().unwrap();
        let k = TIMES
            .iter()
            .position(|x| (x - t).abs() <= 1e-9 * x)
            .unwrap();
        counts.insert(
            (f[0].parse().unwrap(), f[1].to_string(), k),
            f[3].parse().unwrap(),
        );
    }
    Reference {
        molecules: dir.join("g4_events_molecules.csv"),
        edep,
        counts,
    }
}

/// `[irt]` text of one reference event, in Geant4's 2 um water box.
fn event_irt(r: &Reference, event: i64, t_end: f64) -> String {
    let records: Vec<String> = TIMES[1..].iter().map(|t| format!("{t:e}")).collect();
    format!(
        "[irt]\nt_end = {t_end:e}\nstart_contact = false\nbackground_reactions = false\n\
         boundary = \"outflow\"\n\
         record_times = [{}]\n[irt.box]\nlo_nm = -1000\nhi_nm = 1000\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = {T_SEED:e}\nevents = [{event}]\nenergy_ev = {}\n",
        records.join(", "),
        toml_path(&r.molecules),
        r.edep[&event]
    )
}

#[test]
fn irt_single_track_matches_geant4_on_ten_events() {
    let r = reference();
    let dir = scratch("g4");
    let trajectories = 32;
    let mut species = Vec::new();
    // per (species, time): radkin molecules, Geant4 molecules, variance
    let mut sum: HashMap<(usize, usize), (f64, f64, f64)> = HashMap::new();
    let mut init: HashMap<usize, (f64, f64)> = HashMap::new();
    for event in 0..10i64 {
        let cfg = config(&dir, &format!("e{event}.toml"), &event_irt(&r, event, 1e-6));
        let sim = IrtSimulation::from_toml(&cfg).unwrap();
        let out = run_irt_with(&sim, trajectories, 100 + event as u64).unwrap();
        species = out.species.clone();
        // record times: the start, then TIMES[1..]
        assert_eq!(out.record_times.len(), TIMES.len());
        assert_eq!(out.record_times[0], T_SEED);
        for (s, label) in species.iter().enumerate() {
            let g4 = |k: usize| {
                r.counts
                    .get(&(event, label.clone(), k))
                    .copied()
                    .unwrap_or(0.0)
            };
            let e = init.entry(s).or_default();
            e.0 += (g4(0) - out.totals_mean[0][s]).powi(2);
            e.1 += g4(0);
            for k in 1..TIMES.len() {
                let v = sum.entry((s, k)).or_default();
                v.0 += out.totals_mean[k][s];
                v.1 += g4(k);
                v.2 += out.totals_std[k][s].powi(2) * (1.0 + 1.0 / trajectories as f64);
            }
        }
    }
    let energy: f64 = r.edep.values().sum::<f64>() / 100.0;
    let mut worst: f64 = 0.0;
    for (s, label) in species.iter().enumerate() {
        for k in 1..TIMES.len() {
            let (rk, g4, var) = sum[&(s, k)];
            // the two codes start from two realisations of the pre-chemical stage
            let (d2, n0) = init[&s];
            let survive = if n0 > 0.0 { (g4 / n0).min(1.0) } else { 1.0 };
            let sd = (var + d2 * survive * survive).max(1.0).sqrt();
            let z = (rk - g4) / sd;
            worst = worst.max(z.abs());
            assert!(
                z.abs() <= 4.0,
                "{label} at {:e} s: radkin {rk:.1}, Geant4 {g4}, sd {sd:.1} (z {z:+.2})",
                TIMES[k]
            );
        }
    }
    println!("largest |z| on ten events: {worst:.2}");
    // the yields are those of a 100 keV electron track segment
    let g = |label: &str, k: usize| {
        let s = species.iter().position(|x| x == label).unwrap();
        sum[&(s, k)].0 / energy
    };
    assert!((3.8..4.6).contains(&g("e_aq", 2)), "{}", g("e_aq", 2));
    assert!((2.2..3.2).contains(&g("e_aq", 6)), "{}", g("e_aq", 6));
    assert!((0.4..0.9).contains(&g("H2O2", 6)), "{}", g("H2O2", 6));
}

#[test]
fn irt_single_track_same_seed_same_result() {
    let r = reference();
    let dir = scratch("repro");
    let cfg = config(&dir, "e.toml", &event_irt(&r, 3, 1e-7));
    let sim = IrtSimulation::from_toml(&cfg).unwrap();
    let a = run_irt_with(&sim, 6, 7).unwrap();
    let b = run_irt_with(&sim, 6, 7).unwrap();
    assert_eq!(a.totals_mean, b.totals_mean);
    assert_eq!(a.totals_std, b.totals_std);
    assert_eq!(a.stats, b.stats);
    let c = run_irt_with(&sim, 6, 8).unwrap();
    assert_ne!(a.totals_mean, c.totals_mean);
    // a trajectory depends on its stream only, not on the others
    assert_eq!(sim.run_trajectory(7, 2), sim.run_trajectory(7, 2));
    assert_ne!(sim.run_trajectory(7, 2).0, sim.run_trajectory(7, 3).0);
}

#[test]
fn irt_single_track_counts_change_only_through_reactions() {
    let r = reference();
    let dir = scratch("balance");
    // free space: no molecule is lost
    let irt = event_irt(&r, 6, 1e-6).replace("boundary = \"outflow\"\n", "");
    let sim = IrtSimulation::from_toml(config(&dir, "e.toml", &irt)).unwrap();
    let table = SpeciesTable::water();
    let charge: Vec<i64> = sim
        .model
        .species
        .iter()
        .map(|s| i64::from(table.get(s).unwrap().charge))
        .collect();
    for stream in 0..4 {
        let (totals, stats) = sim.run_trajectory(11, stream);
        assert!(stats.pair_reactions > 20, "{stats:?}");
        assert_eq!(stats.absorbed.iter().sum::<u64>(), 0);
        // every pair reaction of the Geant4 network keeps the charge
        let q = |row: &Vec<i64>| -> i64 { row.iter().zip(&charge).map(|(n, c)| n * c).sum() };
        for row in &totals {
            assert_eq!(q(row), q(&totals[0]));
        }
        // the change of every count is the sum of the reactions that fired
        let mut expected = totals[0].clone();
        for p in &sim.model.pairs {
            let n = stats.by_reaction[p.index] as i64;
            expected[p.a] -= n;
            expected[p.b] -= n;
            for &s in &p.products {
                expected[s] += n;
            }
        }
        assert_eq!(&expected, totals.last().unwrap());
        assert_eq!(
            stats.by_reaction.iter().sum::<u64>(),
            stats.pair_reactions + stats.contact_reactions
        );
    }
}

#[test]
fn irt_single_track_no_reaction_keeps_every_count() {
    let dir = scratch("inert");
    // H2 and H2O2 do not react with each other
    let mut csv = String::from("species,x,y,z\n");
    for i in 0..40 {
        let name = if i % 2 == 0 { "H2" } else { "H2O2" };
        csv.push_str(&format!(
            "{name},{},{},0.0\n",
            i as f64 * 0.3,
            (i % 7) as f64
        ));
    }
    let mol = dir.join("inert.csv");
    std::fs::write(&mol, csv).unwrap();
    let irt = format!(
        "[irt]\nt_end = 1e-6\n[[irt.seed]]\nfile = \"{}\"\ntime = 1e-12\n",
        toml_path(&mol)
    );
    let sim = IrtSimulation::from_toml(config(&dir, "inert.toml", &irt)).unwrap();
    let out = run_irt_with(&sim, 3, 1).unwrap();
    assert_eq!(out.stats.pair_reactions + out.stats.contact_reactions, 0);
    assert_eq!(out.record_times.len(), 100);
    let h2 = out.species.iter().position(|s| s == "H2").unwrap();
    let h2o2 = out.species.iter().position(|s| s == "H2O2").unwrap();
    for (mean, std) in out.totals_mean.iter().zip(&out.totals_std) {
        assert_eq!((mean[h2], mean[h2o2]), (20.0, 20.0));
        assert_eq!(mean.iter().sum::<f64>(), 40.0);
        assert!(std.iter().all(|s| *s == 0.0));
    }
}

#[test]
fn irt_single_track_contact_rule_at_the_start() {
    let dir = scratch("contact");
    // H3O+ and OH- 0.1 nm apart: inside the effective radius (1.01 nm) of a
    // fully diffusion-controlled reaction, which leaves no product
    let mol = dir.join("pair.csv");
    std::fs::write(&mol, "species,x,y,z\nH3O_plus,0,0,0\nOH_minus,0.1,0,0\n").unwrap();
    let irt = |contact: bool| {
        format!(
            "[irt]\nt_end = 1e-9\nstart_contact = {contact}\n[[irt.seed]]\nfile = \"{}\"\ntime = 1e-12\n",
            toml_path(&mol)
        )
    };
    let sim = IrtSimulation::from_toml(config(&dir, "on.toml", &irt(true))).unwrap();
    let out = run_irt_with(&sim, 8, 1).unwrap();
    let idx = |l: &str| out.species.iter().position(|s| s == l).unwrap();
    // the record of the start is taken before the reaction, the next one after
    assert_eq!(out.totals_mean[0][idx("H3O_plus")], 1.0);
    assert_eq!(out.totals_mean[1][idx("H3O_plus")], 0.0);
    assert_eq!(out.totals_mean[1][idx("OH_minus")], 0.0);
    assert_eq!(out.stats.contact_reactions, 8);
    assert_eq!(out.stats.pair_reactions, 0);
    // without the rule the pair is looked at again once it has moved apart,
    // and reacts at a sampled time
    let sim = IrtSimulation::from_toml(config(&dir, "off.toml", &irt(false))).unwrap();
    let out = run_irt_with(&sim, 200, 1).unwrap();
    assert_eq!(out.stats.contact_reactions, 0);
    assert!(out.stats.moves >= 200, "{:?}", out.stats);
    assert_eq!(out.totals_mean[1][idx("H3O_plus")], 1.0);
    // After the move the two ions are a few radii apart: they react with the
    // probability r_eff / distance, about one time in two.
    let left = out.totals_mean.last().unwrap()[idx("H3O_plus")];
    assert!((0.3..0.8).contains(&left), "{left}");
    assert_eq!(
        out.stats.pair_reactions,
        (200.0 * (1.0 - left)).round() as u64,
        "{:?}",
        out.stats
    );
}

#[test]
fn irt_single_track_outflow_box_removes_molecules_that_leave() {
    let dir = scratch("outflow");
    // inert molecules in a 20 nm box, moved once at the end: after 1 us they
    // are far outside
    let mut csv = String::from("species,x,y,z\n");
    for i in 0..50 {
        csv.push_str(&format!("H2,{},0,0\n", -9.0 + 0.36 * i as f64));
    }
    let mol = dir.join("h2.csv");
    std::fs::write(&mol, csv).unwrap();
    let irt = |boundary: &str, hi: f64| {
        format!(
            "[irt]\nt_end = 1e-6\nboundary = \"{boundary}\"\n[irt.box]\nlo_nm = -10\nhi_nm = {hi}\n\
             [[irt.seed]]\nfile = \"{}\"\ntime = 1e-12\n",
            toml_path(&mol)
        )
    };
    let h2 = |out: &radkin::irt::IrtOutput| {
        let s = out.species.iter().position(|s| s == "H2").unwrap();
        (out.totals_mean.last().unwrap()[s], out.stats.absorbed[s])
    };
    let sim = IrtSimulation::from_toml(config(&dir, "out.toml", &irt("outflow", 10.0))).unwrap();
    let out = run_irt_with(&sim, 10, 1).unwrap();
    let (left, lost) = h2(&out);
    assert!(left < 1.0, "{left}");
    assert_eq!(lost as f64, 10.0 * (50.0 - left));
    // before the last move nobody has left
    let s = out.species.iter().position(|s| s == "H2").unwrap();
    assert_eq!(out.totals_mean[out.totals_mean.len() - 2][s], 50.0);
    // free space keeps everybody
    let sim = IrtSimulation::from_toml(config(&dir, "free.toml", &irt("none", 10.0))).unwrap();
    assert_eq!(h2(&run_irt_with(&sim, 2, 1).unwrap()), (50.0, 0));
    // a seed outside the box is refused, or dropped with clip
    let e = IrtSimulation::from_toml(config(&dir, "small.toml", &irt("outflow", 5.0)))
        .unwrap_err()
        .to_string();
    assert!(e.contains("outside the box"), "{e}");
    let clipped = irt("outflow", 5.0).replace("[irt]\n", "[irt]\nclip = true\n");
    let sim = IrtSimulation::from_toml(config(&dir, "clip.toml", &clipped)).unwrap();
    assert_eq!(sim.seeds[0].clipped, 11);
    assert_eq!(sim.seeds[0].molecules.len(), 39);
}

#[test]
fn irt_single_track_bad_configurations_are_refused() {
    let dir = scratch("bad");
    let r = reference();
    let error = |name: &str, text: String| {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        IrtSimulation::from_toml(&path).unwrap_err().to_string()
    };
    // a beam and an IRT stage
    let e = error(
        "beam.toml",
        format!(
            "[includes]\ntoml = [\"radkin/config_Labarbe.toml\"]\n{}",
            event_irt(&r, 0, 1e-6)
        ),
    );
    assert!(e.contains("[beam]"), "{e}");
    // no [irt] section
    let e = error(
        "none.toml",
        "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n".into(),
    );
    assert!(e.contains("missing [irt]"), "{e}");
    // a seed at or after the end of the stage
    let late = event_irt(&r, 0, 1e-6)
        .replace("time = 1.001e-12", "time = 2e-6")
        .replace("[irt]\n", "[irt]\nt_start = 1e-12\n");
    let e = error(
        "late.toml",
        format!("[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n{late}"),
    );
    assert!(e.contains("outside the IRT stage"), "{e}");
    // nothing to run
    let e = error(
        "empty.toml",
        "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nt_end = 1e-6\n".into(),
    );
    assert!(e.contains("nothing to run"), "{e}");
    // a species that is not in the chemistry
    let mol = dir.join("x.csv");
    std::fs::write(&mol, "species,x,y,z\nCl^-1,0,0,0\n").unwrap();
    let e = error(
        "species.toml",
        format!(
            "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nt_end = 1e-6\n[[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
            toml_path(&mol)
        ),
    );
    assert!(e.contains("not a dynamic species"), "{e}");
}

#[test]
fn irt_single_track_cli_writes_yields_and_meta() {
    let dir = scratch("cli");
    let out = dir.join("out");
    let run = Command::new(env!("CARGO_BIN_EXE_radkin"))
        .arg("irt")
        .arg(root().join("data/irt_example.toml"))
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(run.status.success(), "{stderr}");
    assert!(stderr.contains("irt run of 66 molecules"), "{stderr}");
    let yields = std::fs::read_to_string(out.join("yields.csv")).unwrap();
    let header = yields.lines().next().unwrap();
    assert!(header.starts_with("t,N_"), "{header}");
    for col in ["N_e_aq", "G_e_aq", "N_e_aq_std", "G_e_aq_std"] {
        assert!(header.split(',').any(|c| c == col), "{col}: {header}");
    }
    assert_eq!(yields.lines().count(), 101);
    // G = N / (E / 100 eV) with the 360 eV of the example
    let cols: Vec<&str> = header.split(',').collect();
    let at = |name: &str| cols.iter().position(|c| *c == name).unwrap();
    let last: Vec<f64> = yields
        .lines()
        .last()
        .unwrap()
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect();
    assert!((last[at("G_e_aq")] - last[at("N_e_aq")] / 3.6).abs() < 1e-9);
    assert!(last[at("N_e_aq")] > 0.0 && last[at("N_e_aq")] < 20.0);
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("irt.json")).unwrap()).unwrap();
    assert_eq!(meta["solver"], "irt");
    assert_eq!(meta["resampling"], "incremental");
    assert_eq!(meta["trajectories"], 32);
    assert_eq!(meta["seed_molecules"], 66);
    assert!(meta["pair_reactions"].as_u64().unwrap() > 100);
    assert_eq!(meta["pair_reaction_types"].as_array().unwrap().len(), 52);
    // an unknown key is an error
    let bad = dir.join("bad.toml");
    std::fs::write(&bad, "[irt]\nt_end = 1e-6\nhandover = 1e-7\n").unwrap();
    let run = Command::new(env!("CARGO_BIN_EXE_radkin"))
        .arg("irt")
        .arg(&bad)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert!(!run.status.success());
    assert!(
        String::from_utf8_lossy(&run.stderr).contains("unknown key 'handover'"),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
}
