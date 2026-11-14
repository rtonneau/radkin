//! Checkpoint and resume of voxel runs (data/voxel_example.toml, shortened).

use std::path::{Path, PathBuf};

use radkin::RadkinError;
use radkin::consts::AVOGADRO;
use radkin::voxel::{Checkpoint, RunControl, VoxelSimulation, write_voxel_dir};

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

/// Coarsening settings that bring the example's track down to level 2 by 3e-9 s
/// (the factor is far above the default so a short run coarsens).
const COARSEN_LEVELS: u32 = 2;
const COARSEN_FACTOR: f64 = 1e5;

/// The example with coarsening on, running to `t_end` with the given snapshots.
fn coarsened(dir: &Path, name: &str, t_end: &str, snapshots: &str) -> VoxelSimulation {
    let mut sim = config(
        dir,
        name,
        SEED,
        &format!("t_end = {t_end}\nsnapshots = {snapshots}\n{RUN}"),
    );
    sim.voxel.max_coarsen_levels = COARSEN_LEVELS;
    sim.voxel.coarsen_tol_factor = COARSEN_FACTOR;
    sim
}

/// Levels of a run's `(t, level)` history.
fn levels(history: &[(f64, u32)]) -> Vec<u32> {
    history.iter().map(|(_, l)| *l).collect()
}

#[test]
fn coarsened_run_resumes_from_a_mid_run_checkpoint() {
    let dir = work_dir("coarse-resume");
    let full = coarsened(&dir, "full.toml", "1e-7", "[1e-8, 3e-8, 1e-7]");
    let reference = full.run().unwrap();
    // Down to level 2 before the checkpoint time, and no change after it.
    assert_eq!(levels(&reference.level_history), [1, 2, 0, 1, 2]);
    assert!(reference.level_history.last().unwrap().0 < 1e-8);

    // "Interrupted" run: its final checkpoint is taken at level 2.
    let mut first = full.clone();
    first.settings.t_end = 1e-8;
    first.settings.snapshots = vec![1e-8];
    let ckpt = dir.join("first.npz");
    let part = first
        .run_with(&RunControl {
            checkpoint: Some(ckpt.clone()),
            resume: None,
        })
        .unwrap();
    assert_eq!(part.level_history, reference.level_history);
    let mid = Checkpoint::read(&ckpt).unwrap();
    assert_eq!(mid.t, 1e-8);

    // The resumed run recomputes the level from the fine-grid state: one level per
    // check, both at the checkpoint time.
    let resumed = full
        .run_with(&RunControl {
            checkpoint: None,
            resume: Some(mid),
        })
        .unwrap();
    assert_eq!(resumed.level_history, [(1e-8, 1), (1e-8, 2)]);
    assert_eq!(resumed.snapshots.len(), 3);
    for (a, b) in resumed.snapshots.iter().zip(&reference.snapshots) {
        assert_eq!(a.t, b.t);
        let d = scaled_diff(&a.state.data, &b.state.data, 1e-6, 1e-12);
        assert!(
            d <= 1.0,
            "t = {:e} s: resumed differs by {d} x tolerance",
            a.t
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// Resume from a checkpoint written at the seed time (damping still pending),
/// with coarsening on or off: matches the uninterrupted run at rtol 1e-6.
fn resume_at_injection_time(coarsen: bool) {
    let dir = work_dir(&format!("inject-resume-{coarsen}"));
    let mut full = coarsened(&dir, "full.toml", "1e-8", "[1e-12, 1e-9, 1e-8]");
    if coarsen {
        // Loose enough that the fresh deposit itself would pass the test.
        full.voxel.coarsen_tol_factor = 1e12;
    } else {
        full.voxel.max_coarsen_levels = 0;
    }
    let reference = full.run().unwrap();
    let mut first = full.clone();
    first.settings.t_end = 1e-12;
    first.settings.snapshots = vec![1e-12];
    let ckpt = dir.join("first.npz");
    first
        .run_with(&RunControl {
            checkpoint: Some(ckpt.clone()),
            resume: None,
        })
        .unwrap();
    let at_seed = Checkpoint::read(&ckpt).unwrap();
    assert_eq!(at_seed.t, 1e-12);
    assert_eq!(at_seed.seeds.len(), 1);
    let resumed = full
        .run_with(&RunControl {
            checkpoint: None,
            resume: Some(at_seed),
        })
        .unwrap();
    // The damped step after the injection runs on the fine grid, as uninterrupted.
    let after: Vec<_> = reference
        .level_history
        .iter()
        .filter(|(t, _)| *t > 1e-12)
        .copied()
        .collect();
    assert_eq!(resumed.level_history, after);
    assert!(
        !coarsen || !after.is_empty(),
        "the reference never coarsened"
    );
    assert_eq!(resumed.snapshots.len(), reference.snapshots.len());
    for (a, b) in resumed.snapshots.iter().zip(&reference.snapshots) {
        assert_eq!(a.t, b.t);
        let d = scaled_diff(&a.state.data, &b.state.data, 1e-6, 1e-12);
        assert!(
            d <= 1.0,
            "t = {:e} s: resumed differs by {d} x tolerance",
            a.t
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn resume_at_injection_time_matches_uninterrupted() {
    resume_at_injection_time(false);
}

#[test]
fn coarsened_resume_at_injection_time_matches_uninterrupted() {
    resume_at_injection_time(true);
}

#[test]
fn snapshots_outputs_and_checkpoints_stay_on_the_fine_grid_when_coarsened() {
    let dir = work_dir("coarse-io");
    let sim = coarsened(&dir, "coarse.toml", "1e-8", "[1e-12, 1e-8]");
    let ckpt = dir.join("coarse.npz");
    let out = sim
        .run_with(&RunControl {
            checkpoint: Some(ckpt.clone()),
            resume: None,
        })
        .unwrap();
    // The last snapshot and the final checkpoint are taken at level 2.
    assert_eq!(levels(&out.level_history), [1, 2, 0, 1, 2]);
    let grid = &sim.voxel.grid;
    let [nx, ny, nz] = grid.dims;
    let n = grid.n_voxels();
    let last = out.snapshots.last().unwrap();
    assert_eq!(last.state.n_voxels, n);
    assert_eq!(last.state.n_species, out.species.len());
    assert_eq!(last.state.data.len(), n * out.species.len());

    // Piecewise constant over the 4x4x4 blocks of level 2.
    let block = 1usize << COARSEN_LEVELS;
    for s in 0..last.state.n_species {
        let f = last.state.field(s);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let head =
                        (i / block * block) + nx * ((j / block * block) + ny * (k / block * block));
                    assert_eq!(
                        f[i + nx * (j + ny * k)],
                        f[head],
                        "species {s} voxel {i},{j},{k}"
                    );
                }
            }
        }
    }

    // The checkpoint holds that same fine-grid state, in the unchanged format.
    let c = Checkpoint::read(&ckpt).unwrap();
    assert_eq!(c.grid, *grid);
    assert_eq!(c.state.n_voxels, n);
    assert_eq!(c.state.data, last.state.data);

    // Output files: fine-grid shapes, the same files as an uncoarsened run.
    let plain = {
        let mut p = sim.clone();
        p.voxel.max_coarsen_levels = 0;
        p.run().unwrap()
    };
    let (d_coarse, d_plain) = (dir.join("o_coarse"), dir.join("o_plain"));
    write_voxel_dir(&d_coarse, &out, grid).unwrap();
    write_voxel_dir(&d_plain, &plain, grid).unwrap();
    let json = |d: &Path| -> serde_json::Value {
        serde_json::from_str(&std::fs::read_to_string(d.join("snapshots.json")).unwrap()).unwrap()
    };
    let (jc, jp) = (json(&d_coarse), json(&d_plain));
    assert_eq!(jc["dims"], serde_json::json!([nx, ny, nz]));
    for key in [
        "dims",
        "spacing_nm",
        "origin_nm",
        "species",
        "unit",
        "layout",
        "snapshots",
    ] {
        assert_eq!(jc[key], jp[key], "{key}");
    }
    let mut names: Vec<_> = std::fs::read_dir(d_coarse.join("snapshots"))
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    names.sort();
    assert_eq!(names.len(), 2 * out.species.len());
    for name in &names {
        let a = std::fs::read(d_coarse.join("snapshots").join(name)).unwrap();
        let b = std::fs::read(d_plain.join("snapshots").join(name)).unwrap();
        assert_eq!(a.len(), b.len(), "{name:?}");
        assert_eq!(
            a[..a.len() - 8 * n],
            b[..b.len() - 8 * n],
            "{name:?} header"
        );
    }
    let csv = |d: &Path| std::fs::read_to_string(d.join("yields.csv")).unwrap();
    let (cc, cp) = (csv(&d_coarse), csv(&d_plain));
    assert_eq!(cc.lines().next(), cp.lines().next());

    // The molecule counts of the last row of yields.csv are those of the fine snapshot.
    let k = 1e-6 * AVOGADRO * grid.voxel_volume_l();
    let row: Vec<f64> = cc
        .lines()
        .last()
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(row[0], last.t);
    for (s, name) in out.species.iter().enumerate() {
        let want = last.state.field(s).iter().sum::<f64>() * k;
        let got = row[1 + s];
        assert!(
            (got - want).abs() <= 1e-9 * want.abs().max(1.0),
            "N_{name}: {got} vs {want}"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn checkpoint_summary_holds_the_coarsening_settings() {
    let dir = work_dir("coarse-summary");
    let sim = coarsened(&dir, "s.toml", "1e-11", "[1e-11]");
    let ckpt = dir.join("s.npz");
    sim.run_with(&RunControl {
        checkpoint: Some(ckpt.clone()),
        resume: None,
    })
    .unwrap();
    let c = Checkpoint::read(&ckpt).unwrap();
    assert_eq!(c.config["max_coarsen_levels"], COARSEN_LEVELS);
    assert_eq!(c.config["coarsen_tol_factor"], COARSEN_FACTOR);

    // A checkpoint written before these keys existed still resumes (the missing
    // keys only show up as a warning), also under other coarsening settings.
    let mut old = c.clone();
    let cfg = old.config.as_object_mut().unwrap();
    cfg.remove("max_coarsen_levels");
    cfg.remove("coarsen_tol_factor");
    let mut longer = sim.clone();
    longer.settings.t_end = 1e-10;
    longer.settings.snapshots = vec![1e-10];
    longer.voxel.max_coarsen_levels = 0;
    let out = longer
        .run_with(&RunControl {
            checkpoint: None,
            resume: Some(old),
        })
        .unwrap();
    assert_eq!(out.snapshots.len(), 1);
    std::fs::remove_dir_all(&dir).ok();
}
