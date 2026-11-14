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
//! ```text
//! cargo run --release --example voxel_bench -- --write-reference
//! cargo run --release --example voxel_bench -- --repeat 3
//! RAYON_NUM_THREADS=1 cargo run --release --example voxel_bench -- --repeat 3
//! ```

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use radkin::voxel::reaction::{reset_stats, stats};
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

fn main() -> radkin::Result<()> {
    let a = Args::parse();
    let mut sim = VoxelSimulation::from_toml(&a.config)?;
    if let Some(t) = a.t_end {
        sim.settings.t_end = t;
    }
    let t_end = sim.settings.t_end;
    sim.settings.snapshots.retain(|t| *t < t_end);
    sim.settings.snapshots.push(t_end);
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

    let mut walls = Vec::new();
    let mut last = None;
    for i in 0..a.repeat.max(1) {
        reset_stats();
        let t0 = Instant::now();
        let out = sim.run()?;
        let wall = t0.elapsed().as_secs_f64();
        let s = stats();
        println!(
            "repeat {i}: {wall:.4} s, {} steps ({} rejected), {} solves",
            out.n_steps, out.n_rejected, s.solves
        );
        walls.push(wall);
        last = Some((out, s));
    }
    let (out, s) = last.expect("at least one repeat");
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
