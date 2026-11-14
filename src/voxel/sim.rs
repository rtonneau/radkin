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
//! checkpoint_every = 2      # checkpoint every 2 snapshot times (optional)
//! ```
//!
//! With a checkpoint path ([`RunControl::checkpoint`]) a checkpoint is
//! written every `checkpoint_every` snapshot times and always at `t_end`.
//! [`RunControl::resume`] starts a run from a checkpoint instead of
//! `t_start`: the grid and species must match, any other setting may differ.

use std::path::{Path, PathBuf};

use crate::chemistry::env::SimEnv;
use crate::config::Config;
use crate::consts::{DEFAULT_ATOL, DEFAULT_RTOL, FIRST_STEP};
use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::solver::{JacobianMode, SolveOptions};
use crate::voxel::checkpoint::{Checkpoint, SeedRecord};
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
    /// Write a checkpoint every this many snapshot times (with a checkpoint
    /// path; the final checkpoint is always written).
    pub checkpoint_every: Option<usize>,
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
        let checkpoint_every = match t.get("checkpoint_every") {
            None => None,
            Some(toml::Value::Integer(i)) if *i >= 1 => Some(*i as usize),
            Some(_) => {
                return Err(RadkinError::Config(
                    "voxel.run.checkpoint_every must be an integer >= 1".into(),
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
            checkpoint_every,
        })
    }

    /// Settings summary stored in checkpoints and compared on resume.
    fn summary(&self) -> serde_json::Value {
        let c = &self.control;
        let dt_max = if c.dt_max.is_finite() {
            serde_json::json!(c.dt_max)
        } else {
            serde_json::Value::Null
        };
        serde_json::json!({
            "t_start": self.t_start,
            "t_end": self.t_end,
            "snapshots": self.snapshots,
            "rtol": c.rtol,
            "atol": c.atol,
            "dt_initial": c.dt_initial,
            "dt_min": c.dt_min,
            "dt_max": dt_max,
            "rannacher_steps": c.rannacher_steps,
            "solver_rtol": self.solve.rtol,
            "solver_atol": self.solve.atol,
            "background_voxels": self.background_voxels,
            "checkpoint_every": self.checkpoint_every,
        })
    }
}

/// Checkpoint options of [`VoxelSimulation::run_with`].
#[derive(Debug, Clone, Default)]
pub struct RunControl {
    /// Where checkpoints are written (overwritten each time).
    pub checkpoint: Option<PathBuf>,
    /// Start from this checkpoint instead of `t_start`.
    pub resume: Option<Checkpoint>,
}

/// FNV-1a hash of `bytes`.
fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(*b)).wrapping_mul(0x0000_0100_0000_01b3)
    })
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
        self.run_with(&RunControl::default())
    }

    /// Configuration summary stored in checkpoints: run settings, a hash of
    /// the chemistry (beam included) and the seeds.
    fn config_summary(&self) -> serde_json::Value {
        let mut v = self.settings.summary();
        v["chemistry_hash"] = format!("{:016x}", fnv(format!("{:?}", self.env).as_bytes())).into();
        v["seeds"] = self
            .seeds
            .iter()
            .map(|s| {
                serde_json::json!({"file": s.spec.file.display().to_string(), "time": s.spec.time})
            })
            .collect();
        v
    }

    /// Seeds still to inject when resuming `ckpt`, and the records of those
    /// already injected.
    fn seeds_after(&self, ckpt: &Checkpoint) -> Result<(Vec<LoadedSeed>, Vec<SeedRecord>)> {
        let mut todo = Vec::new();
        let mut matched = vec![false; ckpt.seeds.len()];
        for s in &self.seeds {
            let rec = SeedRecord {
                time: s.spec.time,
                fingerprint: SeedRecord::fingerprint_of(&s.amount),
            };
            if let Some(i) = (0..ckpt.seeds.len()).find(|i| !matched[*i] && ckpt.seeds[*i] == rec) {
                matched[i] = true;
            } else if s.spec.time < ckpt.t {
                return Err(RadkinError::Checkpoint(format!(
                    "seed {} at {:e} s is before the checkpoint time {:e} s and was not injected",
                    s.spec.file.display(),
                    s.spec.time,
                    ckpt.t
                )));
            } else {
                todo.push(s.clone());
            }
        }
        let n_unmatched = matched.iter().filter(|m| !**m).count();
        if n_unmatched > 0 {
            log::warn!(
                "resume: {n_unmatched} seed(s) injected before the checkpoint are not in the configuration"
            );
        }
        Ok((todo, ckpt.seeds.clone()))
    }

    /// Run with checkpoint writing and/or resuming.
    pub fn run_with(&self, ctl: &RunControl) -> Result<VoxelOutput> {
        let species = self.species();
        let mat = SimMatrices::build(&self.env)?;
        let coeffs: Vec<f64> = species.iter().map(|s| self.voxel.diffusion_of(s)).collect();
        let grid = self.voxel.grid.clone();
        let s = &self.settings;
        let summary = self.config_summary();
        let (state, t0, background, seeds, mut injected) = match &ctl.resume {
            None => {
                let y0 = self.env.initial_values();
                let state = VoxelState::uniform(&y0, grid.n_voxels());
                (state, s.t_start, Some(y0), self.seeds.clone(), Vec::new())
            }
            Some(ckpt) => {
                ckpt.check_compatible(&grid, &species)?;
                if s.t_end <= ckpt.t {
                    return Err(RadkinError::Checkpoint(format!(
                        "t_end ({:e} s) must be after the checkpoint time ({:e} s)",
                        s.t_end, ckpt.t
                    )));
                }
                warn_differences(&ckpt.config, &summary);
                let (todo, done) = self.seeds_after(ckpt)?;
                log::info!(
                    "resuming at t = {:e} s: {} seed(s) already injected, {} to inject",
                    ckpt.t,
                    done.len(),
                    todo.len()
                );
                (
                    ckpt.state.clone(),
                    ckpt.t,
                    ckpt.background.clone(),
                    todo,
                    done,
                )
            }
        };
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(grid.clone(), coeffs),
            state,
            t0,
            s.control.clone(),
            s.solve.clone(),
        );
        if let Some(ckpt) = &ctl.resume {
            run.set_dt(ckpt.dt);
        }
        run.background = background;
        run.skip_background = s.background_voxels;
        run.track_totals();

        // Seeds at or before the current time have been injected (injection
        // comes before the output callback at the same time).
        let mut sorted: Vec<&LoadedSeed> = seeds.iter().collect();
        sorted.sort_by(|a, b| a.spec.time.total_cmp(&b.spec.time));
        let mut n_recorded = 0;
        let mut record_injected = |r: &VoxelRun, injected: &mut Vec<SeedRecord>| {
            while n_recorded < sorted.len() && sorted[n_recorded].spec.time <= r.t {
                let seed = sorted[n_recorded];
                injected.push(SeedRecord {
                    time: seed.spec.time,
                    fingerprint: SeedRecord::fingerprint_of(&seed.amount),
                });
                n_recorded += 1;
            }
        };
        let write_ckpt = |r: &VoxelRun, injected: &[SeedRecord]| -> Result<()> {
            let Some(path) = &ctl.checkpoint else {
                return Ok(());
            };
            Checkpoint {
                t: r.t,
                dt: r.dt(),
                grid: grid.clone(),
                species: species.clone(),
                state: r.state.clone(),
                background: r.background.clone(),
                seeds: injected.to_vec(),
                config: summary.clone(),
            }
            .write(path)?;
            log::info!(
                "checkpoint at t = {:e} s written to {}",
                r.t,
                path.display()
            );
            Ok(())
        };

        let mut snapshots = Vec::new();
        let mut n_outputs = 0usize;
        run_schedule(&mut run, &seeds, &s.snapshots, s.t_end, |r| {
            snapshots.push(Snapshot {
                t: r.t,
                state: r.state.clone(),
            });
            n_outputs += 1;
            if s.checkpoint_every
                .is_some_and(|n| n_outputs.is_multiple_of(n))
            {
                record_injected(r, &mut injected);
                write_ckpt(r, &injected)?;
            }
            Ok(())
        })?;
        record_injected(&run, &mut injected);
        write_ckpt(&run, &injected)?;
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

/// Log every setting that differs between a checkpoint and this run.
fn warn_differences(old: &serde_json::Value, new: &serde_json::Value) {
    let (Some(old), Some(new)) = (old.as_object(), new.as_object()) else {
        return;
    };
    for (k, v) in new {
        match old.get(k) {
            Some(o) if o == v => {}
            o => log::warn!(
                "resume: '{k}' differs from the checkpoint ({} -> {v})",
                o.map_or("absent".to_string(), |o| o.to_string())
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checkpoint_every_setting() {
        let parse = |s: &str| RunSettings::from_table(&toml::from_str(s).unwrap());
        assert_eq!(parse("t_end = 1e-9").unwrap().checkpoint_every, None);
        assert_eq!(
            parse(
                "t_end = 1e-9
checkpoint_every = 3"
            )
            .unwrap()
            .checkpoint_every,
            Some(3)
        );
        assert!(
            parse(
                "t_end = 1e-9
checkpoint_every = 0"
            )
            .is_err()
        );
        assert!(
            parse(
                "t_end = 1e-9
checkpoint_every = 1.5"
            )
            .is_err()
        );
    }

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
