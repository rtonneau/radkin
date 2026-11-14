//! Checkpoint and resume of voxel runs (data/voxel_example.toml, shortened).

use std::path::{Path, PathBuf};

use radkin::RadkinError;
use radkin::voxel::{Checkpoint, RunControl, VoxelSimulation};

/// A fresh directory holding the example's molecule file.
fn work_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin-ckpt-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(
        "data/voxel_example_molecules.csv",
        dir.join("voxel_example_molecules.csv"),
    )
    .unwrap();
    dir
}

/// The example configuration with `[voxel.run]` replaced by `run`, and
/// `seeds` (TOML `[[voxel.seed]]` blocks) in place of its seed.
fn config(dir: &Path, name: &str, seeds: &str, run: &str) -> VoxelSimulation {
    let text = std::fs::read_to_string("data/voxel_example.toml").unwrap();
    let head = &text[..text.find("[[voxel.seed]]").unwrap()];
    let path = dir.join(name);
    std::fs::write(&path, format!("{head}{seeds}\n[voxel.run]\n{run}\n")).unwrap();
    VoxelSimulation::from_toml(&path).unwrap()
}

const SEED: &str = r#"
[[voxel.seed]]
file = "voxel_example_molecules.csv"
time = 1e-12
energy_ev = 360.0
"#;

const RUN: &str = "rtol = 1e-3\natol = 1e-4\ndt_initial = 1e-13";

/// max |a - b| / (atol + rtol |b|) over all values.
fn scaled_diff(a: &[f64], b: &[f64], rtol: f64, atol: f64) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs() / (atol + rtol * y.abs()))
        .fold(0.0, f64::max)
}

#[test]
fn resumed_run_matches_uninterrupted_and_extends() {
    let dir = work_dir("resume");
    let full = config(
        &dir,
        "full.toml",
        SEED,
        &format!("t_end = 1e-9\nsnapshots = [1e-10, 1e-9]\n{RUN}"),
    );
    let ckpt_full = dir.join("full.npz");
    let reference = full
        .run_with(&RunControl {
            checkpoint: Some(ckpt_full.clone()),
            resume: None,
        })
        .unwrap();
    let end = Checkpoint::read(&ckpt_full).unwrap();
    assert_eq!(end.t, 1e-9);
    assert_eq!(end.seeds.len(), 1);
    assert_eq!(
        end.state.data,
        reference.snapshots.last().unwrap().state.data
    );

    // "Interrupted" run: stops at 1e-10 with its final checkpoint.
    let first = config(
        &dir,
        "first.toml",
        SEED,
        &format!("t_end = 1e-10\nsnapshots = [1e-10]\n{RUN}"),
    );
    let ckpt = dir.join("first.npz");
    first
        .run_with(&RunControl {
            checkpoint: Some(ckpt.clone()),
            resume: None,
        })
        .unwrap();
    let mid = Checkpoint::read(&ckpt).unwrap();
    assert_eq!(mid.t, 1e-10);

    // Resume with the full configuration: the seed is not injected twice.
    let resumed = full
        .run_with(&RunControl {
            checkpoint: None,
            resume: Some(mid),
        })
        .unwrap();
    assert_eq!(resumed.snapshots.len(), 2);
    assert_eq!(resumed.snapshots[0].t, 1e-10);
    let a = &resumed.snapshots[1].state.data;
    let b = &reference.snapshots[1].state.data;
    let d = scaled_diff(a, b, 1e-3, 1e-4);
    assert!(d <= 1.0, "resumed run differs from uninterrupted: {d}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn resume_handles_new_and_missing_seeds() {
    let dir = work_dir("seeds");
    let bare = config(&dir, "bare.toml", "", &format!("t_end = 1e-10\n{RUN}"));
    let ckpt = dir.join("bare.npz");
    bare.run_with(&RunControl {
        checkpoint: Some(ckpt.clone()),
        resume: None,
    })
    .unwrap();
    let start = Checkpoint::read(&ckpt).unwrap();
    assert!(start.seeds.is_empty());

    // A seed before the checkpoint that was never injected is an error.
    let late = config(&dir, "late.toml", SEED, &format!("t_end = 1e-9\n{RUN}"));
    let err = late.run_with(&RunControl {
        checkpoint: None,
        resume: Some(start.clone()),
    });
    assert!(matches!(err, Err(RadkinError::Checkpoint(_))), "{err:?}");

    // New conditions: a seed after the checkpoint is injected on top.
    let seed2 = SEED.replace("time = 1e-12", "time = 2e-10");
    let next = config(
        &dir,
        "next.toml",
        &seed2,
        &format!("t_end = 1e-9\nsnapshots = [1e-9]\n{RUN}"),
    );
    let ckpt2 = dir.join("next.npz");
    let out = next
        .run_with(&RunControl {
            checkpoint: Some(ckpt2.clone()),
            resume: Some(start.clone()),
        })
        .unwrap();
    let e_aq = out.species.iter().position(|s| s == "e_aq").unwrap();
    let n = start.state.n_voxels;
    let before: f64 = start.state.data[e_aq * n..(e_aq + 1) * n].iter().sum();
    let after: f64 = out.snapshots[0].state.data[e_aq * n..(e_aq + 1) * n]
        .iter()
        .sum();
    assert!(after > before, "seed not injected: {before} -> {after}");
    assert_eq!(Checkpoint::read(&ckpt2).unwrap().seeds.len(), 1);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn resume_rejects_other_grid_and_past_t_end() {
    let dir = work_dir("grid");
    let sim = config(&dir, "a.toml", "", &format!("t_end = 1e-11\n{RUN}"));
    let ckpt = dir.join("a.npz");
    sim.run_with(&RunControl {
        checkpoint: Some(ckpt.clone()),
        resume: None,
    })
    .unwrap();
    let c = Checkpoint::read(&ckpt).unwrap();

    let err = sim.run_with(&RunControl {
        checkpoint: None,
        resume: Some(c.clone()),
    });
    assert!(matches!(err, Err(RadkinError::Checkpoint(_))), "{err:?}");

    let mut other = sim.clone();
    other.voxel.grid.dims = [8, 8, 10];
    other.settings.t_end = 1e-10;
    let err = other.run_with(&RunControl {
        checkpoint: None,
        resume: Some(c),
    });
    assert!(matches!(err, Err(RadkinError::Checkpoint(_))), "{err:?}");
    std::fs::remove_dir_all(&dir).ok();
}
