//! Multi-grid RDME on an aerated track (ADR 0005): O2, a numerous fast
//! background species, on a coarser level than the radicals.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use radkin::ssa::{RdmeOutput, SsaConfig, run_rdme};
use radkin::voxel::VoxelSimulation;

const K: f64 = 4.0;

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("radkin_ssa_aerated_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Molecule CSV: `n` molecules of each species, uniform in `[lo, hi]^3` nm
/// (deterministic positions).
fn molecules(species: &[&str], n: usize, lo: f64, hi: f64) -> String {
    let mut state = 12345u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut s = String::from("species,x,y,z\n");
    for sp in species {
        for _ in 0..n {
            let p = [0; 3].map(|_| lo + (hi - lo) * next());
            let _ = writeln!(s, "{sp},{},{},{}", p[0], p[1], p[2]);
        }
    }
    s
}

/// Pure water with 250 µmol/L O2 (about 150 O2 per 100 nm voxel), and `n`
/// molecules of each of e_aq, OH, H3O+ and H seeded in `[lo, hi]^3` nm of an
/// 8x8x8 grid of 100 nm voxels.
fn aerated_case(dir: &Path, n: usize, lo: f64, hi: f64) -> PathBuf {
    std::fs::write(
        dir.join("track.csv"),
        molecules(&["e_aq", "OH_r", "H3O_plus", "H_r"], n, lo, hi),
    )
    .unwrap();
    let path = dir.join("run.toml");
    std::fs::write(
        &path,
        r#"
        [includes]
        toml = ["radkin/advanced_water_radiolysis.toml"]

        [concentrations.initial]
        O2 = 250e-6

        [voxel.grid]
        dims = [8, 8, 8]
        spacing_nm = 100.0

        [[voxel.seed]]
        file = "track.csv"
        time = 1e-9

        [voxel.run]
        t_end = 2e-6
        snapshots = [2e-7, 2e-6]
        "#,
    )
    .unwrap();
    path
}

fn ssa(n: usize, o2_level: u32) -> SsaConfig {
    SsaConfig {
        trajectories: n,
        seed: 41,
        species_levels: vec![("O2".into(), o2_level)],
        ..Default::default()
    }
}

/// Largest difference of the species totals over the record times, in
/// combined SEM (floored at the Poisson SEM), and whether every one is
/// within `K` SEM plus one molecule.
fn compare(a: &RdmeOutput, b: &RdmeOutput, n: usize) -> (f64, Option<String>) {
    let mut worst = 0.0f64;
    let mut failure = None;
    for (r, t) in a.record_times.iter().enumerate() {
        for (s, label) in a.species.iter().enumerate() {
            let (ma, mb) = (a.totals_mean[r][s], b.totals_mean[r][s]);
            let var = (a.totals_std[r][s].powi(2) + b.totals_std[r][s].powi(2)) / n as f64;
            let sem = var.max(ma.max(1.0) / n as f64).sqrt();
            worst = worst.max((ma - mb).abs() / sem);
            if (ma - mb).abs() > K * sem + 1.0 && failure.is_none() {
                failure = Some(format!(
                    "t = {t:e} s, {label}: {mb:.1} vs fine {ma:.1} molecules (SEM {sem:.1})"
                ));
            }
        }
    }
    (worst, failure)
}

#[test]
fn ssa_aerated_track_with_coarse_o2_matches_the_fine_grid() {
    // A sparse track: 300 radicals of each kind over 5.4^3 fine voxels, a
    // few per voxel against 150 O2, so O2 stays a background.
    let dir = scratch("track");
    let sim = VoxelSimulation::from_toml(aerated_case(&dir, 300, 130.0, 670.0)).unwrap();
    let n = 16;
    let fine = run_rdme(&sim, &ssa(n, 0)).unwrap();
    let o2 = fine.species.iter().position(|s| s == "O2").unwrap();
    assert!(fine.totals_mean[0][o2] > 50_000.0, "aerated");
    for level in [1, 2] {
        let coarse = run_rdme(&sim, &ssa(n, level)).unwrap();
        assert_eq!(coarse.species_levels[o2], level);
        let (worst, failure) = compare(&fine, &coarse, n);
        println!(
            "aerated track, O2 on level {level}: worst {worst:.2} SEM; hops per trajectory \
             {:.0} vs {:.0} fine ({:.1}x fewer); wall time {:.2} s vs {:.2} s fine",
            coarse.n_hops as f64 / n as f64,
            fine.n_hops as f64 / n as f64,
            fine.n_hops as f64 / coarse.n_hops as f64,
            coarse.exec_time,
            fine.exec_time
        );
        assert!(
            failure.is_none(),
            "O2 on level {level}: {}",
            failure.unwrap()
        );
        assert!(coarse.n_hops * 3 < fine.n_hops);
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// The measurement of `docs/ssa.md` (by hand, about a minute):
/// `cargo test --release --test ssa_rdme_aerated -- --ignored --nocapture`.
#[test]
#[ignore]
fn ssa_aerated_track_measurement() {
    let dir = scratch("measurement");
    let n = 96;
    for (name, count, lo, hi) in [("sparse", 300, 130.0, 670.0), ("dense", 1500, 330.0, 470.0)] {
        let sim = VoxelSimulation::from_toml(aerated_case(&dir, count, lo, hi)).unwrap();
        let fine = run_rdme(&sim, &ssa(n, 0)).unwrap();
        println!(
            "{name}, O2 on level 0: hops per trajectory {:.0}, wall time {:.2} s",
            fine.n_hops as f64 / n as f64,
            fine.exec_time
        );
        for level in [1, 2] {
            let coarse = run_rdme(&sim, &ssa(n, level)).unwrap();
            let (worst, failure) = compare(&fine, &coarse, n);
            println!(
                "{name}, O2 on level {level}: hops per trajectory {:.0} ({:.1}x fewer), wall time \
                 {:.2} s, worst {worst:.2} SEM{}",
                coarse.n_hops as f64 / n as f64,
                fine.n_hops as f64 / coarse.n_hops as f64,
                coarse.exec_time,
                failure
                    .map(|f| format!(" (first failure: {f})"))
                    .unwrap_or_default()
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
