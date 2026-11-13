//! Timing harness used by `tools/benchmark.py`: loads a config and solves it
//! `--repeat` times in-process, then prints the timings as one JSON line.
//!
//! Each repeat builds a fresh `UnitCell` (timed as `load_s`: TOML parsing,
//! includes, matrices) and solves (timed as `solve_s`: `run`, `chunked_run` or
//! `single_pulse_run`). No outputs are written, so file I/O is not measured.
//!
//! ```text
//! cargo run --release --example bench -- CONFIG --t-span 0 1 --mode chunked --repeat 5
//! ```

use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, ValueEnum};
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    Run,
    Chunked,
    SinglePulse,
}

#[derive(Parser)]
struct Args {
    config: PathBuf,
    #[arg(long, num_args = 2, required = true)]
    t_span: Vec<f64>,
    #[arg(long, value_enum, default_value = "run")]
    mode: Mode,
    /// File with output times separated by blanks (mode `run` only).
    #[arg(long)]
    t_eval: Option<PathBuf>,
    #[arg(long)]
    max_step: Option<f64>,
    #[arg(long)]
    max_step_on: Option<f64>,
    #[arg(long)]
    max_step_off: Option<f64>,
    #[arg(long)]
    max_step_trans: Option<f64>,
    #[arg(long)]
    max_step_final: Option<f64>,
    #[arg(long)]
    points_on: Option<usize>,
    #[arg(long)]
    points_off: Option<usize>,
    #[arg(long = "let")]
    let_: Option<f64>,
    #[arg(long, default_value_t = 5)]
    repeat: usize,
}

fn main() -> radkin::Result<()> {
    let a = Args::parse();
    let t_span = (a.t_span[0], a.t_span[1]);
    let t_eval = match &a.t_eval {
        Some(p) => Some(
            std::fs::read_to_string(p)
                .map_err(|e| radkin::RadkinError::io(p, e))?
                .split_whitespace()
                .map(|s| s.parse::<f64>().expect("t_eval value"))
                .collect::<Vec<_>>(),
        ),
        None => None,
    };
    let run_opts = RunOptions {
        max_step: a.max_step,
        t_eval,
        ..Default::default()
    };
    let chunk_opts = ChunkOptions {
        max_step_size_on: a.max_step_on,
        max_step_size_off: a.max_step_off,
        max_step_size_trans: a.max_step_trans,
        max_step_size_final: a.max_step_final,
        points_per_chunk_on: a.points_on,
        points_per_chunk_off: a.points_off,
        ..Default::default()
    };
    let (mut load_s, mut solve_s) = (Vec::new(), Vec::new());
    let (mut nfev, mut nst, mut n_points, mut n_chunks) = (0, 0, 0, 0);
    for _ in 0..a.repeat {
        let t0 = Instant::now();
        let mut cell = UnitCell::from_toml(&a.config)?;
        if let Some(l) = a.let_ {
            cell.change_let(l)?;
        }
        let t1 = Instant::now();
        let res = match a.mode {
            Mode::Run => cell.run(t_span, &run_opts, None)?,
            Mode::Chunked => cell.chunked_run(t_span, &chunk_opts)?,
            Mode::SinglePulse => cell.single_pulse_run(t_span, &chunk_opts)?,
        };
        let t2 = Instant::now();
        load_s.push((t1 - t0).as_secs_f64());
        solve_s.push((t2 - t1).as_secs_f64());
        nfev = res.nfev.iter().sum();
        nst = res.nst.iter().sum();
        n_points = res.t().len();
        n_chunks = res.n_points.len();
    }
    println!(
        "{}",
        serde_json::json!({
            "load_s": load_s,
            "solve_s": solve_s,
            "nfev": nfev,
            "nst": nst,
            "n_points": n_points,
            "n_chunks": n_chunks,
        })
    );
    Ok(())
}
