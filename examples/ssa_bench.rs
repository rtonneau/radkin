//! Exact SSA against tau-leaping on the same case: wall time, events and
//! the agreement of the two ensemble means.
//!
//! ```bash
//! cargo run --release --example ssa_bench -- radicals 1e-12 32
//! cargo run --release --example ssa_bench -- water-const 1e-13 32
//! ```
//!
//! Cases: `radicals` (three reactions, no fast equilibrium), `water-const`
//! (pure water, 1e4 Gy/s for 1 ms) and `water-pulsed` (pure water, 3 pulses
//! of 3 Gy). Arguments: case, volume [L], trajectories, and optionally the
//! tau-leaping epsilon (default 0.03). Set `RAYON_NUM_THREADS=1` for
//! single-core timings.

use radkin::config::Config;
use radkin::gvalues::GValuesBoscolo;
use radkin::physics::Beam;
use radkin::ssa::{EnsembleOptions, SsaConfig, SsaResult, run_ensemble};
use radkin::unit_cell::UnitCell;

const RADICALS: &str = r#"
    [concentrations.initial]
    H_r = 1e-6
    OH_r = 3e-6
    e_aq = 2e-6

    [[reactions.k_reaction]]
    reaction = "H_r -> H2"
    k_value = 3.0e3

    [[reactions.k_reaction]]
    reaction = "2 OH_r -> H2O2"
    k_value = 5.5e9

    [[reactions.k_reaction]]
    reaction = "e_aq + OH_r -> OH_minus"
    k_value = 3.0e10
"#;

fn case(name: &str) -> (UnitCell, (f64, f64)) {
    match name {
        "radicals" => {
            let cfg = Config::from_toml_str(RADICALS, None).unwrap();
            let gvalues = GValuesBoscolo::load(&radkin::consts::data_dir()).unwrap();
            let cell = UnitCell::from_config_with_gvalues(&cfg, gvalues).unwrap();
            (cell, (0.0, 1e-3))
        }
        "water-const" | "water-pulsed" => {
            let mut cell = UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap();
            let (beam, stop) = if name == "water-const" {
                (Beam::constant(1e4, None, None).unwrap(), 1e-3)
            } else {
                (
                    Beam::pulsed_from_dose_per_pulse(3.0, 1e-3, 2e-6, Some(3.0), None, None)
                        .unwrap(),
                    3.5e-3,
                )
            };
            cell.set_beam(beam, "default", Some(2.5)).unwrap();
            (cell, (0.0, stop))
        }
        other => panic!("unknown case '{other}' (radicals, water-const, water-pulsed)"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 && args.len() != 5 {
        eprintln!(
            "usage: ssa_bench <radicals|water-const|water-pulsed> <volume L> <trajectories> [epsilon]"
        );
        std::process::exit(2);
    }
    let volume: f64 = args[2].parse().expect("volume");
    let trajectories: usize = args[3].parse().expect("trajectories");
    let tau_epsilon: f64 = args.get(4).map_or(0.03, |e| e.parse().expect("epsilon"));
    let (cell, t_span) = case(&args[1]);
    let run = |tau_leap: bool| -> SsaResult {
        let cfg = SsaConfig {
            volume: Some(volume),
            trajectories,
            seed: 1,
            tau_leap,
            tau_epsilon,
            rate_correction: None,
            ..Default::default()
        };
        let options = EnsembleOptions {
            t_eval: Some(radkin::unit_cell::thinned_times(t_span.0, t_span.1, 20)),
            ..Default::default()
        };
        run_ensemble(&cell, &cfg, t_span, &options).unwrap()
    };
    let exact = run(false);
    let tau = run(true);
    let n = trajectories as f64;
    println!(
        "case {} | volume {volume:e} L | {trajectories} trajectories | epsilon {tau_epsilon}",
        args[1]
    );
    for (name, r) in [("exact", &exact), ("tau-leap", &tau)] {
        println!(
            "{name:>9}: {:>8.3} s, {:>12.0} events per trajectory, {:>6.1} M events/s",
            r.exec_time,
            r.n_events as f64 / n,
            r.n_events as f64 / r.exec_time / 1e6
        );
    }
    println!("speed-up: {:.1}x", exact.exec_time / tau.exec_time);
    // Largest difference of the two means: in combined standard errors (with
    // a Poisson floor on the variance), and relative to the largest mean
    // of the species over the run.
    let molecule = exact.molecule();
    let n_species = exact.mean[0].len();
    let peak: Vec<f64> = (0..n_species)
        .map(|s| exact.mean.iter().map(|row| row[s]).fold(0.0, f64::max))
        .collect();
    let (mut worst, mut worst_rel) = (0.0f64, 0.0f64);
    for i in 0..exact.t.len() {
        for s in 0..n_species {
            let diff = (exact.mean[i][s] - tau.mean[i][s]).abs();
            let var = (exact.std[i][s].powi(2) + tau.std[i][s].powi(2)) / n;
            let floor = (exact.mean[i][s].max(molecule) * molecule) / n;
            worst = worst.max(diff / var.max(floor).sqrt());
            if peak[s] > 100.0 * molecule {
                worst_rel = worst_rel.max(diff / peak[s]);
            }
        }
    }
    println!(
        "largest difference of the means: {worst:.2} combined SEM, {:.3} % of the species peak",
        100.0 * worst_rel
    );
}
