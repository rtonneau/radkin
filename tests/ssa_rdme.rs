//! Voxel RDME (`radkin voxel` with `solver = "ssa"`) against the
//! deterministic Strang-split voxel run.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use radkin::consts::AVOGADRO;
use radkin::ssa::{RdmeOutput, SsaConfig, run_rdme};
use radkin::voxel::{VoxelOutput, VoxelSimulation};

const K: f64 = 4.0;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_ssa_rdme_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Deterministic pseudo-random numbers in [0, 1) for the molecule positions.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Molecule CSV: `n` molecules of each species, uniform in the box
/// `[lo, hi]^3` nm.
fn molecules(species: &[&str], n: usize, lo: f64, hi: f64) -> String {
    let mut rng = Lcg(12345);
    let mut s = String::from("species,x,y,z\n");
    for sp in species {
        for _ in 0..n {
            let mut p = [0.0; 3];
            for x in &mut p {
                *x = lo + (hi - lo) * rng.next();
            }
            let _ = writeln!(s, "{sp},{},{},{}", p[0], p[1], p[2]);
        }
    }
    s
}

/// Molecules per µmol/L in one voxel.
fn per_umol(sim: &VoxelSimulation) -> f64 {
    1e-6 * AVOGADRO * sim.voxel.grid.voxel_volume_l()
}

fn ssa(trajectories: usize) -> SsaConfig {
    SsaConfig {
        trajectories,
        seed: 31,
        ..Default::default()
    }
}

/// Diffusion only: one species, no reaction.
fn diffusion_case(dir: &Path) -> PathBuf {
    std::fs::write(dir.join("track.csv"), molecules(&["e_aq"], 4000, 5.0, 45.0)).unwrap();
    let path = dir.join("run.toml");
    std::fs::write(
        &path,
        r#"
        [parameters]
        pH = 7

        [reactions.radiolytic]
        e_aq = 2.6

        [voxel.grid]
        dims = [8, 2, 2]
        spacing_nm = 50.0

        [[voxel.seed]]
        file = "track.csv"
        time = 0.0
        energy_ev = 1000.0

        [voxel.run]
        t_end = 5e-6
        snapshots = [2e-7, 1e-6, 5e-6]
        rtol = 1e-5
        atol = 1e-7
        dt_initial = 1e-9
        "#,
    )
    .unwrap();
    path
}

/// Pure water chemistry, radicals seeded in the centre voxel of a 3x3x3 grid
/// of 500 nm voxels (75 molecules per µmol/L and voxel).
fn reaction_case(dir: &Path) -> PathBuf {
    std::fs::write(
        dir.join("track.csv"),
        molecules(&["e_aq", "OH_r", "H3O_plus", "H_r"], 3000, 510.0, 990.0),
    )
    .unwrap();
    let path = dir.join("run.toml");
    std::fs::write(
        &path,
        r#"
        [includes]
        toml = ["radkin/advanced_water_radiolysis.toml"]

        [voxel.grid]
        dims = [3, 3, 3]
        spacing_nm = 500.0

        [[voxel.seed]]
        file = "track.csv"
        time = 1e-9
        energy_ev = 110000.0

        [voxel.run]
        t_end = 2e-4
        snapshots = [1e-6, 1e-5, 5e-5, 2e-4]
        rtol = 1e-5
        atol = 1e-7
        dt_initial = 1e-10
        "#,
    )
    .unwrap();
    path
}

/// Index of `t` in the record times.
fn record(out: &RdmeOutput, t: f64) -> usize {
    out.record_times.iter().position(|x| *x == t).unwrap()
}

/// Domain totals [molecules] of the deterministic run at snapshot `k`.
fn lsoda_totals(sim: &VoxelSimulation, out: &VoxelOutput, k: usize) -> Vec<f64> {
    let unit = per_umol(sim);
    (0..out.species.len())
        .map(|s| out.snapshots[k].state.field(s).iter().sum::<f64>() * unit)
        .collect()
}

#[test]
fn ssa_rdme_diffusion_matches_the_deterministic_field() {
    let dir = scratch("diffusion");
    let sim = VoxelSimulation::from_toml(diffusion_case(&dir)).unwrap();
    let lsoda = sim.run().unwrap();
    let n = 64;
    let out = run_rdme(&sim, &ssa(n)).unwrap();
    assert_eq!(out.n_reactions, 0);
    assert!(out.n_hops > 0);
    let e = out.species.iter().position(|s| s == "e_aq").unwrap();
    let unit = per_umol(&sim);
    let mut worst = 0.0f64;
    for (k, snap) in lsoda.snapshots.iter().enumerate() {
        assert_eq!(out.snapshot_times[k], snap.t);
        // Hops conserve the molecules exactly: no spread on the total.
        let r = record(&out, snap.t);
        assert_eq!(out.totals_mean[r][e], 4000.0);
        assert_eq!(out.totals_std[r][e], 0.0);
        for v in 0..sim.voxel.grid.n_voxels() {
            let ode = snap.state.field(e)[v] * unit;
            let mean = out.mean[k].field(e)[v] * unit;
            let sem = (out.std[k].field(e)[v] * unit / (n as f64).sqrt())
                .max((ode.max(0.0) / n as f64).sqrt());
            // The deterministic field carries its own time-stepping error.
            let allowed = K * sem + 1e-3 * ode + 0.5;
            assert!(
                (mean - ode).abs() <= allowed,
                "t = {:e} s, voxel {v}: rdme {mean:.2} vs lsoda {ode:.2} molecules (SEM {sem:.2})",
                snap.t
            );
            worst = worst.max((mean - ode).abs() / sem);
        }
    }
    println!("diffusion: worst voxel difference {worst:.2} SEM");
    // The field is not uniform at the first snapshot and nearly so at the last.
    let first = lsoda.snapshots[0].state.field(e);
    assert!(first[0] > 50.0 * first[7]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ssa_rdme_reactions_match_the_deterministic_totals() {
    let dir = scratch("reactions");
    let sim = VoxelSimulation::from_toml(reaction_case(&dir)).unwrap();
    let lsoda = sim.run().unwrap();
    let n = 48;
    let out = run_rdme(&sim, &ssa(n)).unwrap();
    assert!(out.n_reactions > 0 && out.n_hops > 0);
    let (mut worst, mut sharp) = (0.0f64, 0);
    for (k, snap) in lsoda.snapshots.iter().enumerate() {
        let r = record(&out, snap.t);
        let ode = lsoda_totals(&sim, &lsoda, k);
        for (s, label) in out.species.iter().enumerate() {
            let mean = out.totals_mean[r][s];
            let sem =
                (out.totals_std[r][s] / (n as f64).sqrt()).max((ode[s].max(0.0) / n as f64).sqrt());
            // One molecule of finite-volume bias, and the tolerance of the
            // deterministic run.
            let allowed = K * sem + 1e-3 * ode[s] + 1.0;
            assert!(
                (mean - ode[s]).abs() <= allowed,
                "t = {:e} s, {label}: rdme {mean:.2} vs lsoda {:.2} molecules (SEM {sem:.2})",
                snap.t,
                ode[s]
            );
            worst = worst.max((mean - ode[s]).abs() / sem);
            if ode[s] > 500.0 {
                sharp += 1;
            }
        }
    }
    println!(
        "reactions: worst total difference {worst:.2} SEM, {sharp} checks above 500 molecules, \
         {} reactions and {} hops per trajectory",
        out.n_reactions / n as u64,
        out.n_hops / n as u64
    );
    assert!(sharp >= 8, "only {sharp} sharp checks");
    // G-values: the net molecules per 100 eV of deposited energy.
    let last = out.record_times.len() - 1;
    let h2o2 = out.species.iter().position(|s| s == "H2O2").unwrap();
    assert!(out.net_mean[last][h2o2] > 100.0);
    assert_eq!(out.seed_energies, vec![(1e-9, 110000.0)]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ssa_rdme_seeds_are_whole_molecules_and_conserved() {
    let dir = scratch("seeds");
    let sim = VoxelSimulation::from_toml(diffusion_case(&dir)).unwrap();
    let out = run_rdme(&sim, &ssa(3)).unwrap();
    let e = out.species.iter().position(|s| s == "e_aq").unwrap();
    // The seed is injected at t_start = 0: 4000 molecules from the first record on.
    assert_eq!(out.record_times[0], 0.0);
    for row in &out.totals_mean {
        assert_eq!(row[e], 4000.0);
    }
    // Net of the initial solution: exactly the seed.
    assert_eq!(out.net_mean[0][e], 4000.0);
    assert_eq!(out.net_std[0][e], 0.0);
    let _ = std::fs::remove_dir_all(&dir);
}

fn radkin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_radkin"))
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

#[test]
fn ssa_rdme_cli_writes_mean_and_std_files() {
    let dir = scratch("cli");
    let cfg = diffusion_case(&dir);
    // Append the solver settings: `[ssa]` after the other tables, and the
    // solver through the flag.
    let mut text = std::fs::read_to_string(&cfg).unwrap();
    text.push_str("\n[ssa]\ntrajectories = 6\nseed = 4\n");
    std::fs::write(&cfg, text).unwrap();
    let run = |out: &Path, threads: &str, extra: &[&str]| -> Output {
        radkin()
            .arg("voxel")
            .arg(&cfg)
            .args(["--solver", "ssa"])
            .args(extra)
            .arg("--out")
            .arg(out)
            .env("RAYON_NUM_THREADS", threads)
            .output()
            .unwrap()
    };
    let (one, many) = (dir.join("one"), dir.join("many"));
    let a = run(&one, "1", &[]);
    assert!(a.status.success(), "{}", stderr(&a));
    let b = run(&many, "4", &[]);
    assert!(b.status.success(), "{}", stderr(&b));

    let yields = std::fs::read_to_string(one.join("yields.csv")).unwrap();
    let header: Vec<&str> = yields.lines().next().unwrap().split(',').collect();
    assert_eq!(header[0], "t");
    for col in ["N_e_aq", "G_e_aq", "N_e_aq_std", "G_e_aq_std"] {
        assert!(header.contains(&col), "{col} missing in {header:?}");
    }
    // The N_ and G_ columns come first, as in the deterministic run.
    let n_species = (header.len() - 1) / 4;
    assert!(
        header[1..=n_species]
            .iter()
            .all(|c| c.starts_with("N_") && !c.ends_with("_std"))
    );
    assert!(
        header[n_species + 1..=2 * n_species]
            .iter()
            .all(|c| c.starts_with("G_"))
    );
    // 4000 e_aq for 1000 eV: G = 400 per 100 eV at every time.
    let g = header.iter().position(|c| *c == "G_e_aq").unwrap();
    for line in yields.lines().skip(1) {
        let v: f64 = line.split(',').nth(g).unwrap().parse().unwrap();
        assert_eq!(v, 400.0);
    }

    for name in ["t000_e_aq.npy", "t000_e_aq_std.npy", "t002_e_aq.npy"] {
        let x = std::fs::read(one.join("snapshots").join(name)).unwrap();
        let y = std::fs::read(many.join("snapshots").join(name)).unwrap();
        assert!(!x.is_empty() && x == y, "{name}: 1 and 4 threads differ");
    }
    assert_eq!(
        yields,
        std::fs::read_to_string(many.join("yields.csv")).unwrap(),
        "1 and 4 threads differ"
    );
    let index: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(one.join("snapshots.json")).unwrap())
            .unwrap();
    assert_eq!(index["solver"], "ssa");
    assert_eq!(index["trajectories"], 6);
    assert_eq!(index["snapshots"].as_array().unwrap().len(), 3);
    assert!(index["n_hops"].as_u64().unwrap() > 0);
    assert!(!one.join("checkpoint.npz").exists());

    // No checkpoint with the ssa solver.
    let c = run(&dir.join("resume"), "1", &["--resume", "nothing.npz"]);
    assert!(!c.status.success());
    assert!(stderr(&c).contains("--resume"), "{}", stderr(&c));
    let _ = std::fs::remove_dir_all(&dir);
}

/// Pure water, radicals seeded in the centre of an 8x8x8 grid of 100 nm
/// voxels, twice; the box homogenises after each seed. `levels` coarsening
/// levels.
fn homogenising_case(dir: &Path, levels: u32) -> PathBuf {
    std::fs::write(
        dir.join("track.csv"),
        molecules(&["e_aq", "OH_r", "H3O_plus"], 1500, 310.0, 490.0),
    )
    .unwrap();
    let path = dir.join(format!("run{levels}.toml"));
    let text = format!(
        r#"
        [includes]
        toml = ["radkin/advanced_water_radiolysis.toml"]

        [voxel]
        max_coarsen_levels = {levels}

        [voxel.grid]
        dims = [8, 8, 8]
        spacing_nm = 100.0

        [[voxel.seed]]
        file = "track.csv"
        time = 1e-9

        [[voxel.seed]]
        file = "track.csv"
        time = 5e-5

        [voxel.run]
        t_end = 1e-4
        snapshots = [1e-6, 1e-5, 1e-4]
        "#
    );
    std::fs::write(&path, text).unwrap();
    path
}

#[test]
fn ssa_rdme_coarsening_matches_the_fine_grid() {
    let dir = scratch("coarsening");
    let fine_sim = VoxelSimulation::from_toml(homogenising_case(&dir, 0)).unwrap();
    let coarse_sim = VoxelSimulation::from_toml(homogenising_case(&dir, 3)).unwrap();
    let n = 24;
    let fine = run_rdme(&fine_sim, &ssa(n)).unwrap();
    let coarse = run_rdme(&coarse_sim, &ssa(n)).unwrap();
    assert!(fine.level_history.is_empty());
    // Trajectory 0: coarsened once the first seed has spread, refined to the
    // fine grid for the second seed, coarsened again afterwards.
    let levels: Vec<u32> = coarse.level_history.iter().map(|(_, l)| *l).collect();
    let refined = levels
        .iter()
        .position(|l| *l == 0)
        .expect("refined for the second seed");
    assert!(
        levels[..refined].iter().any(|l| *l >= 2),
        "levels {levels:?}"
    );
    assert_eq!(coarse.level_history[refined].0, 5e-5);
    assert!(
        levels[refined..].iter().any(|l| *l >= 2),
        "levels {levels:?}"
    );

    let mut worst = 0.0f64;
    for t in &fine.snapshot_times {
        let (rf, rc) = (record(&fine, *t), record(&coarse, *t));
        for (s, label) in fine.species.iter().enumerate() {
            let (mf, mc) = (fine.totals_mean[rf][s], coarse.totals_mean[rc][s]);
            let var =
                (fine.totals_std[rf][s].powi(2) + coarse.totals_std[rc][s].powi(2)) / n as f64;
            let sem = var.max(mf.max(1.0) / n as f64).sqrt();
            assert!(
                (mf - mc).abs() <= K * sem + 1.0,
                "t = {t:e} s, {label}: coarsened {mc:.2} vs fine {mf:.2} molecules (SEM {sem:.2})"
            );
            worst = worst.max((mf - mc).abs() / sem);
        }
    }
    // The snapshots stay on the fine grid and hold every molecule.
    let unit = per_umol(&coarse_sim);
    for (k, t) in coarse.snapshot_times.iter().enumerate() {
        assert_eq!(coarse.mean[k].n_voxels, 512);
        let r = record(&coarse, *t);
        for s in 0..coarse.species.len() {
            let total = coarse.mean[k].field(s).iter().sum::<f64>() * unit;
            assert!(
                (total - coarse.totals_mean[r][s]).abs() <= 1e-9 * total.max(1.0),
                "snapshot {k}, species {s}: {total} vs {}",
                coarse.totals_mean[r][s]
            );
        }
    }
    println!(
        "coarsening: worst total difference {worst:.2} combined SEM; hops per trajectory \
         {} fine, {} coarsened; reactions {} and {}; levels of trajectory 0: {:?}",
        fine.n_hops / n as u64,
        coarse.n_hops / n as u64,
        fine.n_reactions / n as u64,
        coarse.n_reactions / n as u64,
        coarse.level_history
    );
    assert!(coarse.n_hops * 2 < fine.n_hops, "coarsening saves hops");
    let _ = std::fs::remove_dir_all(&dir);
}
