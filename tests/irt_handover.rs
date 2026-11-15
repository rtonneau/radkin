//! `radkin voxel` with an `[irt]` section: the IRT stage handed over to the
//! RDME in the same process.

use std::path::{Path, PathBuf};
use std::process::Command;

use radkin::irt::{IrtSimulation, run_irt_with};
use radkin::ssa::{SolverConfig, handover_seed, run_rdme, run_rdme_from_irt};
use radkin::voxel::VoxelSimulation;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("radkin_irt_handover_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

const HANDOVER: f64 = 1e-7;

/// A voxel run of the Geant4 chemistry on a 2 um box of 32^3 voxels, seeded
/// by the reference event 6, with the `[ssa]` keys `ssa` and the `[irt]`
/// section `irt` (empty: no IRT stage).
fn config(dir: &Path, name: &str, ssa: &str, irt: &str) -> PathBuf {
    let text = format!(
        "solver = \"ssa\"\n[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n\
         [ssa]\ntrajectories = 6\nseed = 4\n{ssa}\n\
         [voxel]\nclip = true\n[voxel.grid]\ndims = [32, 32, 32]\nspacing_nm = 62.5\norigin_nm = [-1000, -1000, -1000]\n\
         [[voxel.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [6]\nenergy_ev = 500.0\n\
         [voxel.run]\nt_end = 1e-6\n{irt}\n",
        toml_path(&root().join("tests/data/irt/g4_events_molecules.csv"))
    );
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

const IRT: &str = "[irt]\nt_end = 1e-7\nboundary = \"reflective\"\nstart_contact = false";

fn load(cfg: &Path) -> (VoxelSimulation, IrtSimulation, radkin::ssa::SsaConfig) {
    (
        VoxelSimulation::from_toml_without_seeds(cfg).unwrap(),
        IrtSimulation::from_toml(cfg).unwrap(),
        SolverConfig::from_toml_file(cfg).unwrap().ssa,
    )
}

#[test]
fn irt_handover_conserves_every_molecule() {
    let dir = scratch("conserve");
    for rule in ["full", "incremental"] {
        let irt_text = format!("{IRT}\nresampling = \"{rule}\"");
        let cfg = config(
            &dir,
            &format!("{rule}.toml"),
            "background_run = false",
            &irt_text,
        );
        let (sim, irt, ssa) = load(&cfg);
        assert!(sim.seeds.is_empty());
        let n_seed = irt.seeds[0].molecules.len();
        assert!(n_seed > 100);
        let out = run_rdme_from_irt(&sim, &irt, &ssa).unwrap();
        let stage = out.irt.as_ref().unwrap();
        assert_eq!(stage.handover_time, HANDOVER);
        assert_eq!(stage.t_start, 1.001e-12);
        assert_eq!(stage.resampling, rule);
        // one time axis: from the start of the stage, through the handover
        assert_eq!(out.record_times[0], 1.001e-12);
        assert_eq!(*out.record_times.last().unwrap(), 1e-6);
        let k_h = out
            .record_times
            .iter()
            .position(|t| *t == HANDOVER)
            .unwrap();
        assert!(k_h > 10 && k_h < out.record_times.len() - 10);

        // the same IRT stage on its own, recorded at the same times
        let mut alone = irt.clone();
        alone.config.record_times = Some(out.record_times[1..k_h].to_vec());
        let reference = run_irt_with(&alone, ssa.trajectories, handover_seed(ssa.seed)).unwrap();
        assert_eq!(reference.record_times, out.record_times[..=k_h].to_vec());
        // before the handover the net molecules are the particles of the stage,
        // and at the handover the grid holds every one of them
        for k in 0..=k_h {
            assert_eq!(
                out.net_mean[k], reference.totals_mean[k],
                "{rule}: record {k} at {:e} s",
                out.record_times[k]
            );
        }
        let particles: f64 = reference.totals_mean[k_h].iter().sum();
        assert_eq!(
            stage.handed_over as f64,
            particles * ssa.trajectories as f64
        );
        assert!(particles > 0.4 * n_seed as f64 && particles < n_seed as f64);
        assert_eq!(stage.stats.pair_reactions, reference.stats.pair_reactions);
        // the totals hold the background of the grid too: 482 H3O+ and OH- at pH 7
        let h3o = out.species.iter().position(|s| s == "H3O_plus").unwrap();
        let background = out.totals_mean[0][h3o] - out.net_mean[0][h3o];
        assert!(
            (background - 481.8).abs() < 4.0 * (481.8f64 / 6.0).sqrt(),
            "{background}"
        );
        // after the handover the grid carries on: reactions fire, molecules
        // are still there, and nothing appears from nowhere
        assert!(out.n_reactions > 0 && out.n_hops > 0);
        let last = out.net_mean.last().unwrap();
        let e = out.species.iter().position(|s| s == "e_aq").unwrap();
        assert!(
            last[e] > 0.0 && last[e] <= out.net_mean[k_h][e],
            "{}",
            last[e]
        );
    }
}

#[test]
fn irt_handover_yields_follow_the_irt_stage_alone() {
    // The RDME after a 100 ns handover against the IRT stage carried on to
    // 1 us: the same chemistry, so the yields stay close.
    let dir = scratch("yields");
    let cfg = config(&dir, "run.toml", "trajectories = 24", IRT);
    let cfg_text = std::fs::read_to_string(&cfg)
        .unwrap()
        .replace("trajectories = 6\n", "");
    std::fs::write(&cfg, cfg_text).unwrap();
    let (sim, irt, ssa) = load(&cfg);
    assert_eq!(ssa.trajectories, 24);
    let out = run_rdme_from_irt(&sim, &irt, &ssa).unwrap();
    let mut alone = irt.clone();
    alone.config.t_end = 1e-6;
    alone.config.record_times = Some(vec![HANDOVER]);
    let reference = run_irt_with(&alone, 96, 9).unwrap();
    let last = out.net_mean.last().unwrap();
    for label in ["e_aq", "OH_r", "H2O2"] {
        let s = out.species.iter().position(|x| x == label).unwrap();
        let (got, want) = (last[s], reference.totals_mean.last().unwrap()[s]);
        assert!(
            (got - want).abs() < 0.25 * want + 3.0,
            "{label}: {got} on the grid against {want} in the IRT stage"
        );
    }
    // G is per 100 eV of the seed
    let csv = radkin::ssa::rdme_yields_csv(&out);
    let header: Vec<&str> = csv.lines().next().unwrap().split(',').collect();
    let col = |name: &str| header.iter().position(|c| *c == name).unwrap();
    let rows: Vec<Vec<&str>> = csv
        .lines()
        .skip(1)
        .map(|l| l.split(',').collect())
        .collect();
    assert_eq!(rows.len(), out.record_times.len());
    let g: f64 = rows[0][col("G_e_aq")].parse().unwrap();
    let e = out.species.iter().position(|x| x == "e_aq").unwrap();
    assert!((g - out.net_mean[0][e] / 5.0).abs() < 1e-9, "{g}");
}

#[test]
fn irt_handover_bad_configurations_are_refused() {
    let dir = scratch("bad");
    let error = |name: &str, ssa: &str, irt: &str| -> String {
        let cfg = config(&dir, name, ssa, irt);
        let sim = VoxelSimulation::from_toml_without_seeds(&cfg).unwrap();
        let ssa = SolverConfig::from_toml_file(&cfg).unwrap().ssa;
        match IrtSimulation::from_toml(&cfg) {
            Err(e) => e.to_string(),
            Ok(irt) => run_rdme_from_irt(&sim, &irt, &ssa).unwrap_err().to_string(),
        }
    };
    // a seed scheduled after the handover time
    let e = error(
        "late.toml",
        "",
        &IRT.replace("t_end = 1e-7", "t_end = 1e-12"),
    );
    assert!(
        e.contains("above the start time") || e.contains("outside the IRT stage"),
        "{e}"
    );
    let late = format!(
        "{IRT}\n[[irt.seed]]\nfile = \"{}\"\ntime = 1.001e-12\nevents = [6]\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 5e-7\nevents = [3]\n",
        toml_path(&root().join("tests/data/irt/g4_events_molecules.csv")),
        toml_path(&root().join("tests/data/irt/g4_events_molecules.csv"))
    );
    let e = error("late2.toml", "", &late);
    assert!(e.contains("outside the IRT stage"), "{e}");
    // walls other than the grid's
    let e = error(
        "free.toml",
        "",
        &IRT.replace("boundary = \"reflective\"\n", ""),
    );
    assert!(e.contains("boundary = \"reflective\""), "{e}");
    let e = error(
        "box.toml",
        "",
        &format!("{IRT}\n[irt.box]\nlo_nm = -1000\nhi_nm = 1200"),
    );
    assert!(e.contains("must be the voxel grid"), "{e}");
    // the ensemble is that of [ssa]
    let e = error("traj.toml", "", &format!("{IRT}\ntrajectories = 3"));
    assert!(e.contains("those of [ssa]"), "{e}");
    // a handover after the end of the run, a snapshot before the handover
    let e = error(
        "after.toml",
        "",
        &IRT.replace("t_end = 1e-7", "t_end = 2e-6"),
    );
    assert!(e.contains("must be before voxel.run.t_end"), "{e}");
    let cfg = config(&dir, "snap.toml", "", IRT);
    let text = std::fs::read_to_string(&cfg).unwrap().replace(
        "[voxel.run]\nt_end = 1e-6",
        "[voxel.run]\nt_end = 1e-6\nsnapshots = [1e-8]",
    );
    std::fs::write(&cfg, text).unwrap();
    let (sim, irt, ssa) = load(&cfg);
    let e = run_rdme_from_irt(&sim, &irt, &ssa).unwrap_err().to_string();
    assert!(e.contains("no field before the handover"), "{e}");
    // seeds left on the grid
    let cfg = config(&dir, "seeds.toml", "", IRT);
    let with_seeds = VoxelSimulation::from_toml(&cfg).unwrap();
    let e = run_rdme_from_irt(&with_seeds, &irt, &ssa)
        .unwrap_err()
        .to_string();
    assert!(e.contains("injected into the IRT stage"), "{e}");
}

#[test]
fn irt_handover_cli() {
    let dir = scratch("cli");
    let radkin = || Command::new(env!("CARGO_BIN_EXE_radkin"));
    let cfg = config(&dir, "run.toml", "", &format!("{IRT}\n[voxel.run2]\n"));
    // `snapshots` inside [voxel.run]: one at the handover, one at the end
    let text = std::fs::read_to_string(&cfg)
        .unwrap()
        .replace("\n[voxel.run2]\n", "")
        .replace(
            "[voxel.run]\nt_end = 1e-6",
            "[voxel.run]\nt_end = 1e-6\nsnapshots = [1e-7, 1e-6]",
        );
    std::fs::write(&cfg, text).unwrap();
    let out = dir.join("out");
    let run = radkin()
        .arg("voxel")
        .arg(&cfg)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(run.status.success(), "{stderr}");
    assert!(
        stderr.contains("irt stage to the handover at 1e-7 s"),
        "{stderr}"
    );
    assert!(
        stderr.contains("ssa voxel run of 32768 voxels, 6 trajectories"),
        "{stderr}"
    );
    let yields = std::fs::read_to_string(out.join("yields.csv")).unwrap();
    let first = yields.lines().nth(1).unwrap();
    assert!(first.starts_with("1.001e-12,"), "{first}");
    assert!(yields.lines().any(|l| l.starts_with("1e-7,")));
    assert!(yields.lines().last().unwrap().starts_with("1e-6,"));
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("snapshots.json")).unwrap())
            .unwrap();
    assert_eq!(meta["irt_stage"]["handover_time"], 1e-7);
    assert_eq!(meta["irt_stage"]["resampling"], "full");
    assert!(meta["irt_stage"]["molecules_handed_over"].as_u64().unwrap() > 300);
    assert_eq!(meta["snapshots"].as_array().unwrap().len(), 2);
    assert!(out.join("snapshots/t000_e_aq.npy").exists());

    // the deterministic solver has no IRT stage
    let run = radkin()
        .arg("voxel")
        .arg(&cfg)
        .args(["--solver", "lsoda"])
        .arg("--out")
        .arg(dir.join("lsoda"))
        .output()
        .unwrap();
    assert!(!run.status.success());
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(stderr.contains("stochastic voxel solver only"), "{stderr}");
}

#[test]
fn irt_handover_absent_section_changes_nothing() {
    // Without [irt] the run is the RDME seeded on the grid, as before.
    let dir = scratch("absent");
    let cfg = config(&dir, "plain.toml", "", "");
    let sim = VoxelSimulation::from_toml(&cfg).unwrap();
    assert_eq!(sim.seeds.len(), 1);
    let ssa = SolverConfig::from_toml_file(&cfg).unwrap().ssa;
    let out = run_rdme(&sim, &ssa).unwrap();
    assert!(out.irt.is_none());
    // the grid starts at voxel.run.t_start = 0 and takes the seed at 1.001 ps
    assert_eq!(out.record_times[0], 0.0);
    let again = run_rdme(&sim, &ssa).unwrap();
    assert_eq!(out.totals_mean, again.totals_mean);
    let files = dir.join("out");
    radkin::ssa::write_rdme_dir(&files, &out, &sim.voxel.grid).unwrap();
    let meta = std::fs::read_to_string(files.join("snapshots.json")).unwrap();
    assert!(!meta.contains("irt_stage"), "{meta}");
}
