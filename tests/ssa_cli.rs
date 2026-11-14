//! `radkin run` with the stochastic solver.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn radkin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_radkin"))
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_ssa_cli_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A config of the Labarbe chemistry (10 Gy at 1e7 Gy/s) with `head` before
/// the tables and `tail` after them.
fn config(dir: &Path, head: &str, tail: &str) -> PathBuf {
    let path = dir.join("run.toml");
    let text = format!("{head}\n[includes]\ntoml = [\"radkin/config_Labarbe.toml\"]\n{tail}\n");
    std::fs::write(&path, text).unwrap();
    path
}

const SSA_TABLE: &str = "[ssa]\nvolume = 1e-16\ntrajectories = 12\nseed = 5";

fn run(cfg: &Path, out: &Path, extra: &[&str], threads: Option<&str>) -> Output {
    let mut cmd = radkin();
    cmd.arg("run")
        .arg(cfg)
        .args(["--t-span", "0", "1e-4"])
        .args(extra)
        .arg("--out")
        .arg(out);
    if let Some(n) = threads {
        cmd.env("RAYON_NUM_THREADS", n);
    }
    cmd.output().unwrap()
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn header(csv: &str) -> Vec<String> {
    csv.lines()
        .next()
        .unwrap()
        .split(',')
        .map(String::from)
        .collect()
}

#[test]
fn ssa_cli_output_does_not_depend_on_the_thread_count() {
    let dir = scratch("threads");
    let cfg = config(&dir, "solver = \"ssa\"", SSA_TABLE);
    let (one, many) = (dir.join("one"), dir.join("many"));
    let a = run(&cfg, &one, &[], Some("1"));
    assert!(a.status.success(), "{}", stderr(&a));
    let b = run(&cfg, &many, &[], Some("4"));
    assert!(b.status.success(), "{}", stderr(&b));
    let csv_one = std::fs::read(one.join("concentrations.csv")).unwrap();
    let csv_many = std::fs::read(many.join("concentrations.csv")).unwrap();
    assert!(csv_one == csv_many, "1 thread and 4 threads differ");

    let meta = |d: &Path| -> serde_json::Value {
        serde_json::from_str(&std::fs::read_to_string(d.join("run.json")).unwrap()).unwrap()
    };
    let (m1, m4) = (meta(&one), meta(&many));
    assert_eq!(m1["solver"], "ssa");
    assert_eq!(m1["volume"], 1e-16);
    assert_eq!(m1["trajectories"], 12);
    assert_eq!(m1["seed"], 5);
    assert_eq!(m1["tau_leap"], false);
    assert!(m1["n_events"].as_u64().unwrap() > 0);
    assert_eq!(m1["n_events"], m4["n_events"]);
    assert!(!one.join("reaction_rates.csv").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ssa_cli_csv_has_the_lsoda_columns_then_the_std_columns() {
    let dir = scratch("layout");
    let cfg = config(&dir, "solver = \"ssa\"", SSA_TABLE);
    let (ssa, lsoda) = (dir.join("ssa"), dir.join("lsoda"));
    let a = run(&cfg, &ssa, &[], None);
    assert!(a.status.success(), "{}", stderr(&a));
    // The flag overrides the config key.
    let b = run(&cfg, &lsoda, &["--solver", "lsoda"], None);
    assert!(b.status.success(), "{}", stderr(&b));
    assert!(lsoda.join("reaction_rates.csv").exists());

    let ssa_csv = std::fs::read_to_string(ssa.join("concentrations.csv")).unwrap();
    let lsoda_csv = std::fs::read_to_string(lsoda.join("concentrations.csv")).unwrap();
    let (hs, hl) = (header(&ssa_csv), header(&lsoda_csv));
    assert_eq!(hs[..hl.len()], hl[..], "LSODA columns first");
    let std_cols = &hs[hl.len()..];
    assert!(!std_cols.is_empty());
    assert!(std_cols.iter().all(|c| c.ends_with("_std")));
    for c in std_cols {
        assert!(hl.contains(&c.trim_end_matches("_std").to_string()), "{c}");
    }
    // 200 log-spaced output times by default, every row complete and finite.
    let rows: Vec<&str> = ssa_csv.lines().skip(1).collect();
    assert_eq!(rows.len(), 200);
    for row in rows {
        let values: Vec<f64> = row.split(',').map(|v| v.parse().unwrap()).collect();
        assert_eq!(values.len(), hs.len());
        assert!(values.iter().all(|v| v.is_finite() && *v >= 0.0));
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ssa_cli_flag_selects_the_solver_and_t_eval_sets_the_times() {
    let dir = scratch("flag");
    // No `solver` key: the flag selects the SSA.
    let cfg = config(&dir, "", SSA_TABLE);
    let times = dir.join("times.txt");
    std::fs::write(&times, "1e-6 1e-5 1e-4").unwrap();
    let out = dir.join("out");
    let r = run(
        &cfg,
        &out,
        &["--solver", "ssa", "--t-eval", times.to_str().unwrap()],
        None,
    );
    assert!(r.status.success(), "{}", stderr(&r));
    let csv = std::fs::read_to_string(out.join("concentrations.csv")).unwrap();
    let t: Vec<f64> = csv
        .lines()
        .skip(1)
        .map(|l| l.split(',').next().unwrap().parse().unwrap())
        .collect();
    assert_eq!(t, [1e-6, 1e-5, 1e-4]);
    // Without the flag the same file runs LSODA.
    let lsoda = dir.join("lsoda");
    let r = run(&cfg, &lsoda, &[], None);
    assert!(r.status.success(), "{}", stderr(&r));
    assert!(lsoda.join("reaction_rates.csv").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ssa_cli_reports_missing_volume_and_lsoda_only_options() {
    let dir = scratch("errors");
    let no_table = config(&dir, "", "");
    let r = run(&no_table, &dir.join("a"), &["--solver", "ssa"], None);
    assert!(!r.status.success());
    assert!(stderr(&r).contains("ssa.volume"), "{}", stderr(&r));

    let cfg = config(&dir, "solver = \"ssa\"", SSA_TABLE);
    let r = run(
        &cfg,
        &dir.join("b"),
        &[
            "--mode",
            "chunked",
            "--rtol",
            "1e-6",
            "--production",
            "H2O2",
        ],
        None,
    );
    assert!(!r.status.success());
    let e = stderr(&r);
    for name in ["--mode", "--rtol", "--production"] {
        assert!(e.contains(name), "{e}");
    }

    let bad = config(&dir, "solver = \"gillespie\"", "");
    let r = run(&bad, &dir.join("c"), &[], None);
    assert!(!r.status.success());
    assert!(stderr(&r).contains("lsoda"), "{}", stderr(&r));
    let _ = std::fs::remove_dir_all(&dir);
}
