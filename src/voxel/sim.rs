//! A complete voxel run from one TOML file: Phase-1 chemistry (with its
//! includes), `[voxel]` grid/diffusion/seeds and `[voxel.run]` settings.
//!
//! ```toml
//! [voxel.run]
//! t_end = 1e-6              # s (required)
//! t_start = 0.0             # s
//! snapshots = [1e-9, 1e-6]  # s
//! rtol = 1e-3               # splitting error control
//! atol = 1e-6               # µmol/L
//! dt_initial = 1e-12
//! dt_max = 1e-7
//! rannacher_steps = 4
//! solver_rtol = 1e-8        # per-voxel rsoda
//! solver_atol = 1e-10       # µmol/L
//! ```

use std::path::Path;

use crate::chemistry::env::SimEnv;
use crate::config::Config;
use crate::consts::{DEFAULT_ATOL, DEFAULT_RTOL, FIRST_STEP};
use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::solver::{JacobianMode, SolveOptions};
use crate::voxel::config::VoxelConfig;
use crate::voxel::diffusion::Diffusion;
use crate::voxel::driver::{StepControl, VoxelRun};
use crate::voxel::grid::VoxelState;
use crate::voxel::seed::{LoadedSeed, load_seeds, run_schedule};

/// `[voxel.run]` settings.
#[derive(Debug, Clone, PartialEq)]
pub struct RunSettings {
    pub t_start: f64,
    pub t_end: f64,
    pub snapshots: Vec<f64>,
    pub control: StepControl,
    pub solve: SolveOptions,
    /// Advance background voxels with one shared integration (default true).
    pub background_voxels: bool,
}

fn opt_num(t: &toml::Table, key: &str) -> Result<Option<f64>> {
    match t.get(key) {
        None => Ok(None),
        Some(toml::Value::Float(f)) => Ok(Some(*f)),
        Some(toml::Value::Integer(i)) => Ok(Some(*i as f64)),
        Some(_) => Err(RadkinError::Config(format!(
            "voxel.run.{key} must be a number"
        ))),
    }
}

impl RunSettings {
    pub fn from_table(t: &toml::Table) -> Result<Self> {
        let t_end = opt_num(t, "t_end")?
            .ok_or_else(|| RadkinError::Config("voxel.run.t_end is required".into()))?;
        let t_start = opt_num(t, "t_start")?.unwrap_or(0.0);
        if t_end.is_nan() || t_end <= t_start {
            return Err(RadkinError::Config(format!(
                "voxel.run: t_end ({t_end}) must be > t_start ({t_start})"
            )));
        }
        let snapshots = match t.get("snapshots") {
            None => Vec::new(),
            Some(v) => v
                .as_array()
                .and_then(|a| {
                    a.iter()
                        .map(|x| x.as_float().or(x.as_integer().map(|i| i as f64)))
                        .collect::<Option<Vec<f64>>>()
                })
                .ok_or_else(|| {
                    RadkinError::Config("voxel.run.snapshots must be a list of numbers".into())
                })?,
        };
        let background_voxels = match t.get("background_voxels") {
            None => true,
            Some(toml::Value::Boolean(b)) => *b,
            Some(_) => {
                return Err(RadkinError::Config(
                    "voxel.run.background_voxels must be a boolean".into(),
                ));
            }
        };
        let d = StepControl::default();
        let control = StepControl {
            rtol: opt_num(t, "rtol")?.unwrap_or(d.rtol),
            atol: opt_num(t, "atol")?.unwrap_or(d.atol),
            dt_initial: opt_num(t, "dt_initial")?.unwrap_or(d.dt_initial),
            dt_min: opt_num(t, "dt_min")?.unwrap_or(d.dt_min),
            dt_max: opt_num(t, "dt_max")?.unwrap_or(d.dt_max),
            rannacher_steps: opt_num(t, "rannacher_steps")?
                .map(|x| x as usize)
                .unwrap_or(d.rannacher_steps),
            ..d
        };
        let solve = SolveOptions {
            rtol: opt_num(t, "solver_rtol")?.unwrap_or(DEFAULT_RTOL),
            atol: opt_num(t, "solver_atol")?.unwrap_or(DEFAULT_ATOL),
            first_step: FIRST_STEP,
            max_step: None,
            use_jac: true,
            jacobian: JacobianMode::Dense,
            t_eval: None,
        };
        Ok(Self {
            t_start,
            t_end,
            snapshots,
            control,
            solve,
            background_voxels,
        })
    }
}

/// One snapshot: the full state at `t`.
#[derive(Debug, Clone)]
pub struct Snapshot {
    pub t: f64,
    pub state: VoxelState,
}

/// Everything a voxel run produces.
#[derive(Debug, Clone)]
pub struct VoxelOutput {
    /// Dynamic species labels, in state order.
    pub species: Vec<String>,
    /// `(t, sums)`: domain sums [µmol/L summed over voxels] per species.
    pub totals: Vec<(f64, Vec<f64>)>,
    /// Background domain sums at the same records as `totals` (empty if
    /// not tracked); G-values are net of them.
    pub background_totals: Vec<Vec<f64>>,
    pub snapshots: Vec<Snapshot>,
    /// `(time, energy_ev)` of every seed with a deposited energy.
    pub seed_energies: Vec<(f64, f64)>,
    pub n_steps: usize,
    pub n_rejected: usize,
}

/// A configured voxel simulation.
#[derive(Debug, Clone)]
pub struct VoxelSimulation {
    pub env: SimEnv,
    pub voxel: VoxelConfig,
    pub settings: RunSettings,
    pub seeds: Vec<LoadedSeed>,
}

impl VoxelSimulation {
    pub fn from_toml(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let cfg = Config::load(path)?;
        let env = SimEnv::from_config(&cfg)?;
        let voxel = VoxelConfig::from_toml_file(path)?;
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        let root: toml::Table =
            toml::from_str(&text).map_err(|e| RadkinError::Config(e.to_string()))?;
        let run = root
            .get("voxel")
            .and_then(|v| v.get("run"))
            .and_then(|v| v.as_table())
            .ok_or_else(|| RadkinError::Config("missing [voxel.run] section".into()))?;
        let settings = RunSettings::from_table(run)?;
        Self::new(env, voxel, settings)
    }

    pub fn new(env: SimEnv, voxel: VoxelConfig, settings: RunSettings) -> Result<Self> {
        for label in voxel.diffusion.keys() {
            if !env.species.contains(label) {
                log::info!("voxel.diffusion: '{label}' is not a species of the chemistry");
            }
        }
        let seeds = load_seeds(&voxel.seeds, &voxel.grid, &env.species, &voxel.policy)?;
        for (i, s) in seeds.iter().enumerate() {
            log::info!(
                "seed {i} at {:e} s: {} molecules placed ({} clipped, {} ignored)",
                s.spec.time,
                s.report.placed,
                s.report.clipped,
                s.report.ignored
            );
        }
        Ok(Self {
            env,
            voxel,
            settings,
            seeds,
        })
    }

    /// Dynamic species labels in state order.
    pub fn species(&self) -> Vec<String> {
        let mut v: Vec<(usize, String)> = self
            .env
            .species
            .dyn_sp
            .values()
            .map(|s| (s.index, s.label.clone()))
            .collect();
        v.sort();
        v.into_iter().map(|(_, l)| l).collect()
    }

    pub fn run(&self) -> Result<VoxelOutput> {
        let species = self.species();
        let mat = SimMatrices::build(&self.env)?;
        let coeffs: Vec<f64> = species.iter().map(|s| self.voxel.diffusion_of(s)).collect();
        let grid = self.voxel.grid.clone();
        let y0 = self.env.initial_values();
        let state = VoxelState::uniform(&y0, grid.n_voxels());
        let s = &self.settings;
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid, coeffs),
            state,
            s.t_start,
            s.control.clone(),
            s.solve.clone(),
        );
        run.background = Some(y0);
        run.skip_background = s.background_voxels;
        run.track_totals();
        let mut snapshots = Vec::new();
        run_schedule(&mut run, &self.seeds, &s.snapshots, s.t_end, |r| {
            snapshots.push(Snapshot {
                t: r.t,
                state: r.state.clone(),
            });
            Ok(())
        })?;
        if run.voxel_steps > 0 {
            log::info!(
                "background voxels: {:.1}% of {} voxel half-steps",
                100.0 * run.background_voxel_steps as f64 / run.voxel_steps as f64,
                run.voxel_steps
            );
        }
        Ok(VoxelOutput {
            species,
            totals: run.totals.take().unwrap_or_default(),
            background_totals: run.background_totals.take().unwrap_or_default(),
            snapshots,
            seed_energies: self
                .seeds
                .iter()
                .filter_map(|s| s.spec.energy_ev.map(|e| (s.spec.time, e)))
                .collect(),
            n_steps: run.steps.len(),
            n_rejected: run.rejected,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn background_voxels_setting() {
        let t: toml::Table = toml::from_str("t_end = 1e-9").unwrap();
        assert!(RunSettings::from_table(&t).unwrap().background_voxels);
        let t: toml::Table = toml::from_str(
            "t_end = 1e-9
background_voxels = false",
        )
        .unwrap();
        assert!(!RunSettings::from_table(&t).unwrap().background_voxels);
        let t: toml::Table = toml::from_str(
            "t_end = 1e-9
background_voxels = 1",
        )
        .unwrap();
        assert!(RunSettings::from_table(&t).is_err());
    }
}
