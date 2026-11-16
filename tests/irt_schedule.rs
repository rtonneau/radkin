//! Tracks injected into the IRT stage at scheduled times.

use std::path::{Path, PathBuf};

use radkin::irt::{IrtSimulation, irt_yields_csv, run_irt_with};

/// Allowed distance in standard errors.
const K: f64 = 4.0;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_irt_sched_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

fn load(dir: &Path, name: &str, text: &str) -> IrtSimulation {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    IrtSimulation::from_toml(&path).unwrap()
}

/// `[[irt.seed]]` of the reference event 3 at `time`, moved by `offset_um`
/// along x.
fn track(time: f64, offset_um: f64, energy_ev: f64) -> String {
    format!(
        "[[irt.seed]]\nfile = \"{}\"\ntime = {time:e}\nevents = [3]\noffset_nm = [{}, 0, 0]\nenergy_ev = {energy_ev}\n",
        toml_path(&root().join("tests/data/irt/g4_events_molecules.csv")),
        offset_um * 1000.0
    )
}

const HEAD: &str = "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nstart_contact = false\nbackground_reactions = false\n";

#[test]
fn irt_schedule_widely_spaced_tracks_reproduce_the_single_track() {
    let dir = scratch("spaced");
    let t0 = 1.001e-12;
    let gap = 2e-7;
    let t_end = 1e-6;
    let trajectories = 96;
    // one track: its counts 1 us, 0.8 us and 0.6 us after it started
    let single = load(
        &dir,
        "single.toml",
        &format!(
            "{HEAD}t_end = {:e}\nrecord_times = [{:e}, {:e}]\n{}",
            t0 + t_end,
            t0 + t_end - gap,
            t0 + t_end - 2.0 * gap,
            track(t0, 0.0, 400.0)
        ),
    );
    let one = run_irt_with(&single, trajectories, 61).unwrap();
    assert_eq!(one.record_times.len(), 4);
    // three copies of it, 50 um apart, 0.2 us one after the other
    let three = load(
        &dir,
        "three.toml",
        &format!(
            "{HEAD}t_end = {:e}\nrecord_times = [{:e}]\n{}{}{}",
            t0 + t_end,
            t0 + 0.5 * gap,
            track(t0, 0.0, 400.0),
            track(t0 + gap, 50.0, 400.0),
            track(t0 + 2.0 * gap, 100.0, 400.0)
        ),
    );
    let n_molecules = single.seeds[0].molecules.len();
    assert!(n_molecules > 50);
    assert_eq!(three.seeds.len(), 3);
    let all = run_irt_with(&three, trajectories, 62).unwrap();
    // records: the start, half a gap later, the end (seed times are not records)
    assert_eq!(all.record_times.len(), 3);
    let sem = |o: &radkin::irt::IrtOutput, k: usize, s: usize| {
        o.totals_std[k][s] / (trajectories as f64).sqrt()
    };
    let last = all.record_times.len() - 1;
    let mut checked = 0;
    for s in 0..all.species.len() {
        // the single track at its ages 0.6, 0.8 and 1 us: records 1, 2, 3
        let expected: f64 = (1..=3).map(|k| one.totals_mean[k][s]).sum();
        let var: f64 =
            (1..=3).map(|k| sem(&one, k, s).powi(2)).sum::<f64>() + sem(&all, last, s).powi(2);
        let got = all.totals_mean[last][s];
        assert!(
            (got - expected).abs() <= K * var.sqrt() + 1e-9,
            "{}: {got} against {expected} (sd {})",
            all.species[s],
            var.sqrt()
        );
        if expected > 5.0 {
            checked += 1;
        }
    }
    assert!(checked >= 4, "{checked}");
    // the three tracks started with three times the molecules of one
    let injected: i64 = three.seeds.iter().map(|s| s.molecules.len() as i64).sum();
    assert_eq!(injected, 3 * n_molecules as i64);
}

#[test]
fn irt_schedule_injection_times_offsets_and_energies() {
    let dir = scratch("times");
    let text = format!(
        "{HEAD}t_end = 1e-6\nt_start = 1e-9\nrecord_times = [5e-9, 2e-8, 5e-8]\n{}{}",
        track(1e-8, 0.0, 400.0),
        track(3e-8, 30.0, 250.0)
    );
    let sim = load(&dir, "two.toml", &text);
    let n = sim.seeds[0].molecules.len() as f64;
    // the offset moves every molecule of the second seed
    for (a, b) in sim.seeds[0].molecules.iter().zip(&sim.seeds[1].molecules) {
        assert_eq!(a.0, b.0);
        assert!((b.1[0] - a.1[0] - 30e-6).abs() < 1e-15);
        assert_eq!((a.1[1], a.1[2]), (b.1[1], b.1[2]));
    }
    let out = run_irt_with(&sim, 4, 63).unwrap();
    // the start, 5 ns, the first seed, 20 ns, the second seed, 50 ns, the end
    assert_eq!(
        out.record_times,
        vec![1e-9, 5e-9, 2e-8, 5e-8, 1e-6],
        "seed times are not record times"
    );
    let total = |k: usize| out.totals_mean[k].iter().sum::<f64>();
    assert_eq!(total(0), 0.0);
    assert_eq!(total(1), 0.0);
    // one track at 20 ns (some molecules have reacted), two at 50 ns
    assert!(total(2) > 0.5 * n && total(2) <= n, "{}", total(2));
    assert!(total(3) > total(2) + 0.5 * n, "{} {}", total(3), total(2));
    assert_eq!(out.seed_energies, vec![(1e-8, 400.0), (3e-8, 250.0)]);
    // G divides by the energy injected so far
    let csv = irt_yields_csv(&out);
    let header: Vec<&str> = csv.lines().next().unwrap().split(',').collect();
    let col = |name: &str| header.iter().position(|c| *c == name).unwrap();
    let rows: Vec<Vec<&str>> = csv
        .lines()
        .skip(1)
        .map(|l| l.split(',').collect())
        .collect();
    let value = |k: usize, name: &str| rows[k][col(name)].parse::<f64>().unwrap();
    assert_eq!(rows[0][col("G_e_aq")], "nan");
    assert!((value(2, "G_e_aq") - value(2, "N_e_aq") / 4.0).abs() < 1e-12);
    assert!((value(3, "G_e_aq") - value(3, "N_e_aq") / 6.5).abs() < 1e-12);
    assert!((value(4, "G_e_aq_std") - value(4, "N_e_aq_std") / 6.5).abs() < 1e-12);
}

#[test]
fn irt_schedule_injected_molecules_join_the_pair_search() {
    let dir = scratch("join");
    // an e_aq at the start and an OH_r 50 ns later, in a box of 8 nm with
    // reflective walls: they meet within a few tens of ns
    let e = dir.join("e.csv");
    std::fs::write(&e, "species,x,y,z\ne_aq,2,2,2\n").unwrap();
    let oh = dir.join("oh.csv");
    std::fs::write(&oh, "species,x,y,z\nOH_r,0,0,0\n").unwrap();
    let text = format!(
        "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = 2.95e10\nirt_type = \"full\"\n\
         [irt]\nt_end = 3e-7\nboundary = \"reflective\"\nstart_contact = false\nrecord_times = [4e-8]\n\
         [irt.box]\nlo_nm = 0\nhi_nm = 8\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 0\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 5e-8\noffset_nm = [6, 6, 6]\n",
        toml_path(&e),
        toml_path(&oh)
    );
    let sim = load(&dir, "join.toml", &text);
    let out = run_irt_with(&sim, 200, 64).unwrap();
    let idx = |l: &str| out.species.iter().position(|s| s == l).unwrap();
    // before the OH_r arrives the e_aq is alone
    assert_eq!(out.totals_mean[1][idx("e_aq")], 1.0);
    assert_eq!(out.totals_mean[1][idx("OH_r")], 0.0);
    // 250 ns later nearly every pair has reacted (rate about 1e8 s-1)
    let last = out.totals_mean.last().unwrap();
    assert!(last[idx("e_aq")] < 0.05, "{}", last[idx("e_aq")]);
    assert_eq!(last[idx("e_aq")], last[idx("OH_r")]);
    assert_eq!(last[idx("OH_minus")], 1.0 - last[idx("e_aq")]);
    // a seed after the end of the stage is an error
    let late = text.replace("time = 5e-8", "time = 4e-7");
    let path = dir.join("late.toml");
    std::fs::write(&path, late).unwrap();
    let err = IrtSimulation::from_toml(&path).unwrap_err().to_string();
    assert!(err.contains("outside the IRT stage"), "{err}");
}
