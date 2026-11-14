//! Benchmark and accuracy gate of the voxel reaction-diffusion run.
//!
//! Runs a voxel case (default: `data/voxel_example.toml`, seeded by
//! `data/voxel_example_molecules.csv`) `--repeat` times in-process and prints the
//! median wall time, the splitting steps, the rsoda solve counts and the accuracy
//! of the final state against a tight-tolerance reference.
//!
//! The reference is the same case with the splitting error control tightened to
//! `--reference-rtol` / `--reference-atol` (defaults 1e-6, 1e-6 umol/L). It is
//! written once by `--write-reference` and read by every normal run; the file
//! stores the case (grid, species, `t_end`) and the final per-voxel state. A
//! normal run compares its own final state with it and prints:
//!
//! - the relative error of every species' domain total (sum over voxels);
//! - the max per-voxel scaled error `|y - y_ref| / (atol + rtol * |y_ref|)` with
//!   the run's own control `rtol` / `atol` (below 1 means within the tolerance
//!   the run asked for).
//!
//! The last stdout line is a stable `key=value` summary starting with
//! `VOXEL_BENCH`, for scripts.
//!
//! With `--max-coarsen-levels N` (N > 0, overriding `[voxel] max_coarsen_levels`) the run
//! is coarsened, and the same case is also run on the fine grid (`--repeat` times, same
//! settings otherwise) for comparison. The bench then prints both times, the level
//! history, the relative error of every species' domain total against the fine run and
//! the max per-voxel difference to the fine run scaled by the merge tolerance
//! (`coarsen_tol_factor * (atol + rtol * |y_fine|)`, below 1 = inside what a merge
//! allows), and adds `fine_median_s`, `speedup`, `coarse_levels`, `coarse_max_total_err`
//! and `coarse_max_scaled_err` to the last line. `--coarsen-tol-factor F` overrides
//! `[voxel] coarsen_tol_factor`. Use a late `--t-end` (e.g. 1e-5) so that the case
//! homogenises.
//!
//! ```text
//! cargo run --release --example voxel_bench -- --write-reference
//! cargo run --release --example voxel_bench -- --repeat 3
//! RAYON_NUM_THREADS=1 cargo run --release --example voxel_bench -- --repeat 3
//! cargo run --release --example voxel_bench -- --t-end 1e-5 --max-coarsen-levels 2
//! ```

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use radkin::voxel::reaction::{ReactionStats, reset_stats, stats};
use radkin::voxel::{VoxelOutput, VoxelSimulation, VoxelState};
use serde::{Deserialize, Serialize};

const DEFAULT_REFERENCE: &str =
    ".scratch/tests/2026-11-14__reduce-computational-burden/voxel_reference.json";

#[derive(Parser)]
struct Args {
    /// Voxel TOML configuration.
    #[arg(default_value = "data/voxel_example.toml")]
    config: PathBuf,
    /// Number of timed runs.
    #[arg(long, default_value_t = 3)]
    repeat: usize,
    /// End time [s] (overrides `[voxel.run] t_end`, 1e-7 s in the default case).
    #[arg(long)]
    t_end: Option<f64>,
    /// Reference file (read for the comparison, written by `--write-reference`).
    #[arg(long, default_value = DEFAULT_REFERENCE)]
    reference: PathBuf,
    /// Run the tight-tolerance reference, write it to `--reference` and exit.
    #[arg(long)]
    write_reference: bool,
    /// Splitting `rtol` of the reference run.
    #[arg(long, default_value_t = 1e-6)]
    reference_rtol: f64,
    /// Splitting `atol` of the reference run [umol/L].
    #[arg(long, default_value_t = 1e-6)]
    reference_atol: f64,
    /// Coarsening levels (overrides `[voxel] max_coarsen_levels`). Above 0 the
    /// fine grid is also run and compared.
    #[arg(long)]
    max_coarsen_levels: Option<u32>,
    /// Tolerance factor of the coarsening test (overrides `[voxel] coarsen_tol_factor`).
    #[arg(long)]
    coarsen_tol_factor: Option<f64>,
}

#[derive(Serialize, Deserialize)]
struct Reference {
    species: Vec<String>,
    dims: [usize; 3],
    t_end: f64,
    rtol: f64,
    atol: f64,
    n_steps: usize,
    wall_s: f64,
    /// Final state, species-major (`data[s * n_voxels + v]`) [umol/L].
    data: Vec<f64>,
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Final snapshot of a run (the bench always snapshots at `t_end`).
fn final_state(out: &VoxelOutput, t_end: f64) -> &VoxelState {
    let s = out.snapshots.last().expect("no snapshot");
    assert!(
        (s.t - t_end).abs() <= 1e-12 * t_end,
        "last snapshot at {} s, not t_end",
        s.t
    );
    &s.state
}

/// Wall times and last output and counters of one case.
type Timed = (Vec<f64>, VoxelOutput, ReactionStats);

/// Run every case of `sims` `repeat` times, alternating between them within each repeat
/// (so that warm-up and clock drift hit all cases alike), printing each wall time.
fn timed_runs(sims: &[(&str, &VoxelSimulation)], repeat: usize) -> radkin::Result<Vec<Timed>> {
    let mut res: Vec<Timed> = Vec::new();
    for i in 0..repeat.max(1) {
        for (k, (label, sim)) in sims.iter().enumerate() {
            reset_stats();
            let t0 = Instant::now();
            let out = sim.run()?;
            let wall = t0.elapsed().as_secs_f64();
            let s = stats();
            println!(
                "{label} {i}: {wall:.4} s, {} steps ({} rejected), {} solves",
                out.n_steps, out.n_rejected, s.solves
            );
            if i == 0 {
                res.push((vec![wall], out, s));
            } else {
                res[k].0.push(wall);
                res[k].1 = out;
                res[k].2 = s;
            }
        }
    }
    Ok(res)
}

fn main() -> radkin::Result<()> {
    let a = Args::parse();
    let mut sim = VoxelSimulation::from_toml(&a.config)?;
    if let Some(t) = a.t_end {
        sim.settings.t_end = t;
    }
    let t_end = sim.settings.t_end;
    sim.settings.snapshots.retain(|t| *t < t_end);
    sim.settings.snapshots.push(t_end);
    if let Some(l) = a.max_coarsen_levels {
        sim.voxel.max_coarsen_levels = l;
    }
    if let Some(f) = a.coarsen_tol_factor {
        sim.voxel.coarsen_tol_factor = f;
    }
    let dims = sim.voxel.grid.dims;
    let n_vox = sim.voxel.grid.n_voxels();
    let threads = std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));

    if a.write_reference {
        let mut r = sim.clone();
        r.settings.control.rtol = a.reference_rtol;
        r.settings.control.atol = a.reference_atol;
        let t0 = Instant::now();
        let out = r.run()?;
        let wall_s = t0.elapsed().as_secs_f64();
        let reference = Reference {
            species: out.species.clone(),
            dims,
            t_end,
            rtol: a.reference_rtol,
            atol: a.reference_atol,
            n_steps: out.n_steps,
            wall_s,
            data: final_state(&out, t_end).data.clone(),
        };
        if let Some(dir) = a.reference.parent() {
            std::fs::create_dir_all(dir).map_err(|e| radkin::RadkinError::io(dir, e))?;
        }
        let text = serde_json::to_string(&reference).expect("serialize reference");
        std::fs::write(&a.reference, text).map_err(|e| radkin::RadkinError::io(&a.reference, e))?;
        println!(
            "reference: rtol {:e}, atol {:e}, {} steps ({} rejected), {wall_s:.3} s, written to {}",
            a.reference_rtol,
            a.reference_atol,
            out.n_steps,
            out.n_rejected,
            a.reference.display()
        );
        return Ok(());
    }

    // With coarsening, the fine grid runs alongside (alternating) for the comparison.
    let mut fine_sim = sim.clone();
    fine_sim.voxel.max_coarsen_levels = 0;
    let mut timed = if sim.voxel.max_coarsen_levels > 0 {
        timed_runs(&[("repeat", &sim), ("fine repeat", &fine_sim)], a.repeat)?
    } else {
        timed_runs(&[("repeat", &sim)], a.repeat)?
    };
    let fine_timed = (timed.len() > 1).then(|| timed.remove(1));
    let (walls, out, s) = timed.remove(0);
    let med = median(&mut walls.clone());
    let min = walls.iter().copied().fold(f64::INFINITY, f64::min);
    let c = &sim.settings.control;
    println!(
        "case {}: grid {}x{}x{} ({n_vox} voxels), t_end {t_end:e} s, control rtol {:e} atol {:e}, solver rtol {:e} atol {:e}, {threads} threads",
        a.config.display(),
        dims[0],
        dims[1],
        dims[2],
        c.rtol,
        c.atol,
        sim.settings.solve.rtol,
        sim.settings.solve.atol
    );
    println!(
        "time: median {med:.4} s, min {min:.4} s over {} repeats",
        walls.len()
    );
    println!(
        "steps: {} accepted, {} rejected; voxel half-steps {} ({} background); rsoda solves {}; nfev {}, nst {}",
        out.n_steps,
        out.n_rejected,
        s.voxel_steps,
        s.background_voxel_steps,
        s.solves,
        s.nfev,
        s.nst
    );

    // Final species totals (sum over voxels) and spread, to record a baseline and to
    // check that a case homogenises: spread = (max - min) / mean over voxels.
    {
        let st = final_state(&out, t_end);
        println!("final species totals (sum over voxels [umol/L]) and spread (max-min)/mean:");
        let mut max_spread = 0.0_f64;
        let mut line = String::new();
        for (s_i, name) in out.species.iter().enumerate() {
            let v = &st.data[s_i * n_vox..(s_i + 1) * n_vox];
            let sum: f64 = v.iter().sum();
            let (lo, hi) = v
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                    (l.min(*x), h.max(*x))
                });
            let mean = sum / n_vox as f64;
            let spread = if mean.abs() > 0.0 {
                (hi - lo) / mean.abs()
            } else {
                0.0
            };
            if sum.abs() > 1e-12 {
                max_spread = max_spread.max(spread);
            }
            println!("  {name:>10}: total {sum:.9e}, spread {spread:.3e}");
            line.push_str(&format!(" total_{name}={sum:.9e}"));
        }
        println!("max spread over species with a nonzero total: {max_spread:.3e}");
        println!("VOXEL_TOTALS{line}");
    }

    let mut acc = String::from(" acc=none");
    match std::fs::read_to_string(&a.reference) {
        Err(e) => println!(
            "accuracy: no reference at {} ({e}); run with --write-reference",
            a.reference.display()
        ),
        Ok(text) => {
            let r: Reference = serde_json::from_str(&text).expect("parse reference");
            if r.species != out.species || r.dims != dims || r.t_end != t_end {
                println!(
                    "accuracy: reference does not match the case (species, grid or t_end); rewrite it with --write-reference"
                );
            } else {
                let st = final_state(&out, t_end);
                let mut worst_total = (0.0_f64, "");
                let mut totals = Vec::new();
                for (s_i, name) in out.species.iter().enumerate() {
                    let run = &st.data[s_i * n_vox..(s_i + 1) * n_vox];
                    let rf = &r.data[s_i * n_vox..(s_i + 1) * n_vox];
                    let (sr, sf): (f64, f64) = (run.iter().sum(), rf.iter().sum());
                    // Floor the denominator at atol so empty species do not blow up.
                    let err = (sr - sf).abs() / sf.abs().max(c.atol);
                    if err > worst_total.0 {
                        worst_total = (err, name);
                    }
                    totals.push((name.as_str(), err));
                }
                let max_scaled = st
                    .data
                    .iter()
                    .zip(&r.data)
                    .map(|(y, yr)| (y - yr).abs() / (c.atol + c.rtol * yr.abs()))
                    .fold(0.0, f64::max);
                println!(
                    "accuracy vs reference (rtol {:e}, atol {:e}, {} steps, {:.3} s):",
                    r.rtol, r.atol, r.n_steps, r.wall_s
                );
                for (name, err) in &totals {
                    println!("  total rel. error {name:>10}: {err:.3e}");
                }
                println!(
                    "max domain-total rel. error: {:.3e} ({}); max per-voxel scaled error: {max_scaled:.3e}",
                    worst_total.0, worst_total.1
                );
                acc = format!(
                    " max_total_err={:.6e} max_scaled_err={max_scaled:.6e}",
                    worst_total.0
                );
            }
        }
    }
    if let Some(fine) = &fine_timed {
        acc.push_str(&coarse_comparison(&sim, &out, &walls, s.voxel_steps, fine));
    }
    println!(
        "VOXEL_BENCH median_s={med:.6} min_s={min:.6} repeats={} threads={threads} voxels={n_vox} steps={} rejected={} voxel_steps={} background_voxel_steps={} solves={} nfev={} nst={}{acc}",
        walls.len(),
        out.n_steps,
        out.n_rejected,
        s.voxel_steps,
        s.background_voxel_steps,
        s.solves,
        s.nfev,
        s.nst
    );
    Ok(())
}

/// Compare the coarsened `out` with the fine run `fine` (see the module doc).
/// Returns the `key=value` fields to add to the summary line.
fn coarse_comparison(
    sim: &VoxelSimulation,
    out: &VoxelOutput,
    walls: &[f64],
    coarse_voxel_steps: u64,
    fine: &Timed,
) -> String {
    let t_end = sim.settings.t_end;
    let n_vox = sim.voxel.grid.n_voxels();
    let c = &sim.settings.control;
    let (fine_walls, fine_out, fs) = fine;
    let med = median(&mut walls.to_vec());
    let fine_med = median(&mut fine_walls.clone());
    let fine_min = fine_walls.iter().copied().fold(f64::INFINITY, f64::min);
    let min = walls.iter().copied().fold(f64::INFINITY, f64::min);
    println!(
        "coarsening: max_coarsen_levels {}, coarsen_tol_factor {:e}",
        sim.voxel.max_coarsen_levels, sim.voxel.coarsen_tol_factor
    );
    println!(
        "time: fine median {fine_med:.4} s (min {fine_min:.4}), coarsened median {med:.4} s (min {min:.4}), speedup x{:.3} (min/min x{:.3})",
        fine_med / med,
        fine_min / min
    );
    println!(
        "fine: {} steps, voxel half-steps {}; coarsened: {} steps, voxel half-steps {coarse_voxel_steps} ({:.1}% of fine)",
        fine_out.n_steps,
        fs.voxel_steps,
        out.n_steps,
        100.0 * coarse_voxel_steps as f64 / fs.voxel_steps as f64
    );
    let hist: Vec<String> = out
        .level_history
        .iter()
        .map(|(t, l)| format!("{t:.3e}:{l}"))
        .collect();
    println!(
        "level history ({} changes): [{}]",
        out.level_history.len(),
        hist.join(", ")
    );
    let st = final_state(out, t_end);
    let fst = final_state(fine_out, t_end);
    let mut worst_total = (0.0_f64, "");
    for (s_i, name) in out.species.iter().enumerate() {
        let range = s_i * n_vox..(s_i + 1) * n_vox;
        let (sr, sf): (f64, f64) = (
            st.data[range.clone()].iter().sum(),
            fst.data[range].iter().sum(),
        );
        let err = (sr - sf).abs() / sf.abs().max(c.atol);
        println!("  total rel. error vs fine {name:>10}: {err:.3e}");
        if err > worst_total.0 {
            worst_total = (err, name);
        }
    }
    let scaled = |scale: f64| {
        st.data
            .iter()
            .zip(&fst.data)
            .map(|(y, yf)| (y - yf).abs() / (scale * (c.atol + c.rtol * yf.abs())))
            .fold(0.0, f64::max)
    };
    let max_scaled = scaled(sim.voxel.coarsen_tol_factor);
    let max_scaled_control = scaled(1.0);
    println!(
        "vs fine run: max domain-total rel. error {:.3e} ({}); max per-voxel difference {max_scaled:.3e} merge tolerances ({max_scaled_control:.3e} control tolerances)",
        worst_total.0, worst_total.1
    );
    format!(
        " fine_median_s={fine_med:.6} speedup={:.4} coarse_levels={} coarse_max_total_err={:.6e} coarse_max_scaled_err={max_scaled:.6e}",
        fine_med / med,
        out.level_history.len(),
        worst_total.0
    )
}
