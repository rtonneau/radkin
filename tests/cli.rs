use std::process::Command;

fn radkin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_radkin"))
}

#[test]
fn config_without_beam_fails_with_a_message() {
    let out = radkin()
        .args([
            "run",
            "data/advanced_water_radiolysis.toml",
            "--t-span",
            "1e-7",
            "1e-3",
            "--out",
        ])
        .arg(std::env::temp_dir().join("radkin_cli_nobeam"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("No beam"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn labarbe_run_writes_every_file() {
    let dir = std::env::temp_dir().join(format!("radkin_cli_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let out = radkin()
        .args([
            "run",
            "data/config_Labarbe.toml",
            "--t-span",
            "1e-7",
            "1e-3",
        ])
        .args(["--production", "H2O2", "--consumption", "O2", "--out"])
        .arg(&dir)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    for f in [
        "concentrations.csv",
        "reaction_rates.csv",
        "species_integrals.csv",
        "production_H2O2.csv",
        "consumption_O2.csv",
        "run.json",
    ] {
        assert!(dir.join(f).exists(), "{f} missing");
    }
    let conc = std::fs::read_to_string(dir.join("concentrations.csv")).unwrap();
    let header: Vec<&str> = conc.lines().next().unwrap().split(',').collect();
    assert_eq!(header[0], "t");
    assert_eq!(header.len(), 22 + 7 + 1);
    let rates = std::fs::read_to_string(dir.join("reaction_rates.csv")).unwrap();
    assert!(
        rates.lines().next().unwrap().ends_with("EZ_85"),
        "{}",
        rates.lines().next().unwrap()
    );
    let meta: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("run.json")).unwrap()).unwrap();
    assert!(meta["g_values"].as_object().unwrap().contains_key("G_e_aq"));
    assert_eq!(
        meta["reaction_strings"]["G_e_aq"],
        "RadiolyticReaction('e_aq')"
    );
}

/// Last row of `concentrations.csv` of a Labarbe run with extra `flags`.
fn labarbe_final_row(tag: &str, flags: &[&str]) -> Vec<f64> {
    let dir = std::env::temp_dir().join(format!("radkin_cli_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let out = radkin()
        .args([
            "run",
            "data/config_Labarbe.toml",
            "--t-span",
            "1e-7",
            "1e-3",
        ])
        .args(flags)
        .arg("--out")
        .arg(&dir)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let conc = std::fs::read_to_string(dir.join("concentrations.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    conc.lines()
        .last()
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect()
}

#[test]
fn sparse_jacobian_flag_matches_dense_run() {
    let dense = labarbe_final_row("dense", &[]);
    let sparse = labarbe_final_row("sparse", &["--sparse-jacobian"]);
    assert_eq!(dense.len(), sparse.len());
    let peak = dense.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (i, (s, d)) in sparse.iter().zip(&dense).enumerate() {
        assert!(
            (s - d).abs() <= 1e-6 * d.abs() + 1e-10 * peak,
            "column {i}: {s} vs {d}"
        );
    }
}

#[test]
fn voxel_resume_from_checkpoint() {
    let dir = std::env::temp_dir().join(format!("radkin_cli_resume_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let first = dir.join("first");
    let out = radkin()
        .args(["voxel", "data/voxel_example.toml", "--out"])
        .arg(&first)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ckpt = first.join("checkpoint.npz");
    assert!(ckpt.exists());
    // Same directory for the resumed run is refused.
    let out = radkin()
        .args(["voxel", "data/voxel_example.toml", "--resume"])
        .arg(&ckpt)
        .arg("--out")
        .arg(&first)
        .output()
        .unwrap();
    assert!(!out.status.success());
    // The checkpoint is at t_end: resuming the same config has nothing to do.
    let second = dir.join("second");
    let out = radkin()
        .args(["voxel", "data/voxel_example.toml", "--resume"])
        .arg(&ckpt)
        .arg("--out")
        .arg(&second)
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("t_end"));
    // Extension to a later t_end.
    std::fs::copy(
        "data/voxel_example_molecules.csv",
        dir.join("voxel_example_molecules.csv"),
    )
    .unwrap();
    let text = std::fs::read_to_string("data/voxel_example.toml")
        .unwrap()
        .replace("t_end = 1e-7", "t_end = 2e-7");
    std::fs::write(dir.join("longer.toml"), text).unwrap();
    let out = radkin()
        .arg("voxel")
        .arg(dir.join("longer.toml"))
        .arg("--resume")
        .arg(&ckpt)
        .arg("--out")
        .arg(&second)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(second.join("checkpoint.npz").exists());
    let yields = std::fs::read_to_string(second.join("yields.csv")).unwrap();
    let last_t: f64 = yields
        .lines()
        .last()
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!((last_t - 2e-7).abs() < 1e-18, "{last_t}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn voxel_example_writes_yields_and_snapshots() {
    let dir = std::env::temp_dir().join(format!("radkin_cli_voxel_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let out = radkin()
        .args(["voxel", "data/voxel_example.toml", "--out"])
        .arg(&dir)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let yields = std::fs::read_to_string(dir.join("yields.csv")).unwrap();
    let header = yields.lines().next().unwrap();
    assert!(header.starts_with("t,") && header.contains("N_e_aq") && header.contains("G_OH_r"));
    let index: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("snapshots.json")).unwrap())
            .unwrap();
    let snaps = index["snapshots"].as_array().unwrap();
    assert_eq!(snaps.len(), 4);
    let first = snaps[0]["files"][0].as_str().unwrap();
    let bytes = std::fs::read(dir.join(first)).unwrap();
    assert_eq!(&bytes[..6], b"\x93NUMPY");
    let _ = std::fs::remove_dir_all(&dir);
}
