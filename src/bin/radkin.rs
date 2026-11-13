//! `radkin` command line: run a configuration and write CSV/JSON outputs.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand, ValueEnum};
use radkin::RadkinError;
use radkin::output::{RunMeta, write_run_dir};
use radkin::solver::JacobianMode;
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell};

#[derive(Parser)]
#[command(
    version,
    about = "Homogeneous water radiolysis kinetics (Rust port of radiopyo)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    /// One integration, the beam read at every evaluation (`UnitCell.run`).
    Run,
    /// Chunks cut at the beam transitions (`prepare_chunked_run`).
    Chunked,
    /// Single-pulse chunks (`prepare_single_pulse`).
    SinglePulse,
}

#[derive(Subcommand)]
enum Command {
    /// Run a simulation and write its outputs.
    Run(Box<RunArgs>),
    /// Voxel reaction-diffusion run seeded by Geant4-DNA molecules
    /// (`[voxel]` sections of the configuration).
    Voxel(VoxelArgs),
}

#[derive(clap::Args)]
struct VoxelArgs {
    /// TOML configuration file with `[voxel.grid]`, `[[voxel.seed]]` and `[voxel.run]`.
    config: PathBuf,
    /// Output directory.
    #[arg(long)]
    out: PathBuf,
}

fn voxel(a: &VoxelArgs) -> radkin::Result<()> {
    let sim = radkin::voxel::VoxelSimulation::from_toml(&a.config)?;
    let started = std::time::Instant::now();
    let out = sim.run()?;
    radkin::voxel::write_voxel_dir(&a.out, &out, &sim.voxel.grid)?;
    eprintln!(
        "radkin: voxel run of {} voxels, {} steps ({} rejected) in {:.1} s; {} snapshots written to {}",
        sim.voxel.grid.n_voxels(),
        out.n_steps,
        out.n_rejected,
        started.elapsed().as_secs_f64(),
        out.snapshots.len(),
        a.out.display()
    );
    Ok(())
}

#[derive(clap::Args)]
struct RunArgs {
    /// TOML configuration file.
    config: PathBuf,
    /// Start and end times [s].
    #[arg(long, num_args = 2, value_names = ["START", "STOP"], required = true)]
    t_span: Vec<f64>,
    #[arg(long, value_enum, default_value = "run")]
    mode: Mode,
    /// File with output times separated by blanks (mode `run` only).
    #[arg(long)]
    t_eval: Option<PathBuf>,
    #[arg(long)]
    rtol: Option<f64>,
    #[arg(long)]
    atol: Option<f64>,
    #[arg(long)]
    max_step: Option<f64>,
    #[arg(long)]
    first_step: Option<f64>,
    #[arg(long)]
    max_step_on: Option<f64>,
    #[arg(long)]
    max_step_off: Option<f64>,
    #[arg(long)]
    max_step_trans: Option<f64>,
    #[arg(long)]
    max_step_final: Option<f64>,
    #[arg(long)]
    trans_time: Option<f64>,
    #[arg(long)]
    points_on: Option<usize>,
    #[arg(long)]
    points_off: Option<usize>,
    /// Label of the config beam to use.
    #[arg(long)]
    beam: Option<String>,
    /// LET [keV/µm] (rewrites the Boscolo G-values).
    #[arg(long = "let")]
    let_: Option<f64>,
    /// Temperature [degC].
    #[arg(long)]
    temperature: Option<f64>,
    /// Finite-difference Jacobian instead of the analytic one.
    #[arg(long)]
    no_jac: bool,
    /// Sparse Jacobian and sparse LU instead of dense ones.
    #[arg(long)]
    sparse_jacobian: bool,
    /// Species whose gross production is written (comma separated).
    #[arg(long, value_delimiter = ',')]
    production: Vec<String>,
    /// Species whose gross consumption is written (comma separated).
    #[arg(long, value_delimiter = ',')]
    consumption: Vec<String>,
    /// Output directory.
    #[arg(long)]
    out: PathBuf,
}

fn read_times(path: &PathBuf) -> radkin::Result<Vec<f64>> {
    let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
    text.split_whitespace()
        .map(str::parse::<f64>)
        .collect::<Result<_, _>>()
        .map_err(|e| RadkinError::Value(format!("{}: {e}", path.display())))
}

fn run(a: &RunArgs) -> radkin::Result<()> {
    let mut cell = UnitCell::from_toml(&a.config)?;
    if let Some(t) = a.temperature {
        cell.set_temperature(t)?;
    }
    if let Some(label) = &a.beam {
        cell.use_beam(label, a.let_)?;
    } else if let Some(l) = a.let_ {
        cell.beam()?;
        cell.change_let(l)?;
    }
    let t_span = (a.t_span[0], a.t_span[1]);
    let t_eval = a.t_eval.as_ref().map(read_times).transpose()?;
    if !matches!(a.mode, Mode::Run) && t_eval.is_some() {
        return Err(RadkinError::Value(
            "--t-eval is only available with --mode run".into(),
        ));
    }
    let run_opts = RunOptions {
        atol: a.atol,
        rtol: a.rtol,
        max_step: a.max_step,
        first_step: a.first_step,
        use_jac: !a.no_jac,
        jacobian: if a.sparse_jacobian {
            JacobianMode::Sparse
        } else {
            JacobianMode::Dense
        },
        t_eval,
        ..Default::default()
    };
    let chunk_opts = ChunkOptions {
        run: RunOptions {
            t_eval: None,
            max_step: None,
            ..run_opts.clone()
        },
        max_step_size_on: a.max_step_on,
        max_step_size_off: a.max_step_off,
        max_step_size_trans: a.max_step_trans,
        max_step_size_final: a.max_step_final,
        trans_time: a.trans_time,
        points_per_chunk_on: a.points_on,
        points_per_chunk_off: a.points_off,
        y0: None,
    };
    let (res, mode) = match a.mode {
        Mode::Run => (cell.run(t_span, &run_opts, None)?, "run"),
        Mode::Chunked => (cell.chunked_run(t_span, &chunk_opts)?, "chunked"),
        Mode::SinglePulse => (cell.single_pulse_run(t_span, &chunk_opts)?, "single-pulse"),
    };
    let config = a.config.display().to_string();
    let meta = RunMeta::from_result(&res, &config, mode, t_span, cell.beams.current_name());
    write_run_dir(&a.out, &res, &meta, &a.production, &a.consumption)?;
    eprintln!(
        "radkin: {} points written to {}",
        res.t().len(),
        a.out.display()
    );
    Ok(())
}

fn main() -> ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let cli = Cli::parse();
    let result = match &cli.command {
        Command::Run(a) => run(a),
        Command::Voxel(a) => voxel(a),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("radkin: error: {e}");
            ExitCode::FAILURE
        }
    }
}
