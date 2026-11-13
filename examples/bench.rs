//! Timing harness used by `tools/benchmark.py`: loads a config and solves it
//! `--repeat` times in-process, then prints the timings as one JSON line.
//!
//! Each repeat builds a fresh `UnitCell` (timed as `load_s`: TOML parsing,
//! includes, matrices) and solves (timed as `solve_s`: `run`, `chunked_run` or
//! `single_pulse_run`). No outputs are written, so file I/O is not measured.
//!
//! `--jacobian sparse` solves with the sparse Jacobian instead of the dense
//! one; `--jacobian both` times both and adds a `"sparse"` object (its
//! `solve_s`, `nfev`, `nst`) to the dense line.
//!
//! `--breakdown` (needs the `timing` feature) adds a `"breakdown"` object to
//! each line: per repeat, the time in kernel RHS and Jacobian calls and in
//! rsoda's linear-solver `prepare` and `solve` calls, with their counts (see
//! `radkin::solver::timing`).
//!
//! ```text
//! cargo run --release --example bench -- CONFIG --t-span 0 1 --mode chunked --repeat 5
//! cargo run --release --example bench -- CONFIG --t-span 0 1 --jacobian both
//! cargo run --release --features timing --example bench -- CONFIG --t-span 0 1 --breakdown
//! ```

use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, ValueEnum};
use radkin::result::OdeResult;
use radkin::solver::JacobianMode;
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    Run,
    Chunked,
    SinglePulse,
}

#[derive(Clone, Copy, PartialEq, ValueEnum)]
enum Jacobian {
    Dense,
    Sparse,
    Both,
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
    #[arg(long, value_enum, default_value = "dense")]
    jacobian: Jacobian,
    /// Report the time breakdown of each solve (needs the `timing` feature).
    #[arg(long)]
    breakdown: bool,
}

struct Timings {
    load_s: Vec<f64>,
    solve_s: Vec<f64>,
    breakdown: Vec<serde_json::Value>,
    last: Option<OdeResult>,
}

/// The breakdown of the solves since the last call, as JSON (seconds).
#[cfg(feature = "timing")]
fn take_breakdown() -> serde_json::Value {
    let b = radkin::solver::timing::take();
    serde_json::json!({
        "n_rhs": b.n_rhs, "t_rhs": b.t_rhs.as_secs_f64(),
        "n_jac": b.n_jac, "t_jac": b.t_jac.as_secs_f64(),
        "n_prepare": b.n_prepare, "t_prepare": b.t_prepare.as_secs_f64(),
        "n_solve": b.n_solve, "t_solve": b.t_solve.as_secs_f64(),
    })
}

#[cfg(not(feature = "timing"))]
fn take_breakdown() -> serde_json::Value {
    serde_json::Value::Null
}

fn time_runs(
    a: &Args,
    t_span: (f64, f64),
    run_opts: &RunOptions,
    chunk_opts: &ChunkOptions,
) -> radkin::Result<Timings> {
    let mut t = Timings {
        load_s: Vec::new(),
        solve_s: Vec::new(),
        breakdown: Vec::new(),
        last: None,
    };
    for _ in 0..a.repeat {
        let t0 = Instant::now();
        let mut cell = UnitCell::from_toml(&a.config)?;
        if let Some(l) = a.let_ {
            cell.change_let(l)?;
        }
        take_breakdown();
        let t1 = Instant::now();
        let res = match a.mode {
            Mode::Run => cell.run(t_span, run_opts, None)?,
            Mode::Chunked => cell.chunked_run(t_span, chunk_opts)?,
            Mode::SinglePulse => cell.single_pulse_run(t_span, chunk_opts)?,
        };
        let t2 = Instant::now();
        t.load_s.push((t1 - t0).as_secs_f64());
        t.solve_s.push((t2 - t1).as_secs_f64());
        if a.breakdown {
            t.breakdown.push(take_breakdown());
        }
        t.last = Some(res);
    }
    Ok(t)
}

fn main() -> radkin::Result<()> {
    let a = Args::parse();
    if a.breakdown && cfg!(not(feature = "timing")) {
        return Err(radkin::RadkinError::Value(
            "--breakdown needs the `timing` feature (cargo run --features timing ...)".into(),
        ));
    }
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
    let first = if a.jacobian == Jacobian::Sparse {
        JacobianMode::Sparse
    } else {
        JacobianMode::Dense
    };
    let options = |jacobian: JacobianMode| {
        let run_opts = RunOptions {
            max_step: a.max_step,
            t_eval: t_eval.clone(),
            jacobian,
            ..Default::default()
        };
        let chunk_opts = ChunkOptions {
            run: RunOptions {
                jacobian,
                ..Default::default()
            },
            max_step_size_on: a.max_step_on,
            max_step_size_off: a.max_step_off,
            max_step_size_trans: a.max_step_trans,
            max_step_size_final: a.max_step_final,
            points_per_chunk_on: a.points_on,
            points_per_chunk_off: a.points_off,
            ..Default::default()
        };
        (run_opts, chunk_opts)
    };
    let counts = |r: &Option<OdeResult>| match r {
        Some(r) => (
            r.nfev.iter().sum::<usize>(),
            r.nst.iter().sum::<usize>(),
            r.t().len(),
            r.n_points.len(),
        ),
        None => (0, 0, 0, 0),
    };
    let (ro, co) = options(first);
    let main = time_runs(&a, t_span, &ro, &co)?;
    let (nfev, nst, n_points, n_chunks) = counts(&main.last);
    let mut line = serde_json::json!({
        "load_s": main.load_s,
        "solve_s": main.solve_s,
        "nfev": nfev,
        "nst": nst,
        "n_points": n_points,
        "n_chunks": n_chunks,
    });
    if a.breakdown {
        line["breakdown"] = main.breakdown.into();
    }
    if a.jacobian == Jacobian::Both {
        let (ro, co) = options(JacobianMode::Sparse);
        let sparse = time_runs(&a, t_span, &ro, &co)?;
        let (nfev, nst, _, _) = counts(&sparse.last);
        line["sparse"] = serde_json::json!({ "solve_s": sparse.solve_s, "nfev": nfev, "nst": nst });
        if a.breakdown {
            line["sparse"]["breakdown"] = sparse.breakdown.into();
        }
    }
    println!("{line}");
    Ok(())
}
