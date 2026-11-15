//! A complete IRT run from one TOML file: the chemistry (with its
//! includes), the `[irt]` section and its seeds.
//!
//! Every trajectory starts from the same seed molecules with its own random
//! stream (config seed, stream index = trajectory index), as the SSA
//! ensemble does.

use std::path::Path;
use std::time::Instant;

use crate::chemistry::env::SimEnv;
use crate::config::Config;
use crate::error::{RadkinError, Result};
use crate::irt::config::{BoundaryKind, IrtConfig};
use crate::irt::model::IrtModel;
use crate::irt::sample::Vec3;
use crate::irt::solver::{Boundary, Stage, StageOptions, StageStats};
use crate::matrices::SimMatrices;
use crate::ssa::ensemble::{Moments, map_ordered};
use crate::ssa::rng::SsaRng;
use crate::unit_cell::thinned_times;
use crate::voxel::seed::{SeedPolicy, SeedSpec, read_molecules};

/// Number of record times when `irt.record_times` is not given.
pub const N_RECORDS: usize = 100;

/// A seed ready to inject: its molecules as particles.
#[derive(Debug, Clone, PartialEq)]
pub struct IrtSeed {
    pub spec: SeedSpec,
    /// `(species, position [m])` of every molecule.
    pub molecules: Vec<(u32, Vec3)>,
    /// Molecules of the file dropped by the ignore list.
    pub ignored: usize,
    /// Molecules outside the box, dropped (`clip`).
    pub clipped: usize,
}

/// Read a seed file into particles of the dynamic species of `species`.
pub fn load_seed(
    spec: &SeedSpec,
    species: &[String],
    policy: &SeedPolicy,
    domain: Option<(Vec3, Vec3)>,
) -> Result<IrtSeed> {
    let mut mols = read_molecules(&spec.file)?;
    if let Some(ev) = &spec.events {
        mols.retain(|m| m.event_id.is_some_and(|e| ev.contains(&e)));
    }
    let mut seed = IrtSeed {
        spec: spec.clone(),
        molecules: Vec::with_capacity(mols.len()),
        ignored: 0,
        clipped: 0,
    };
    for m in &mols {
        if policy.ignore.contains(&m.species) {
            seed.ignored += 1;
            continue;
        }
        let label = policy
            .aliases
            .get(&m.species)
            .map(String::as_str)
            .unwrap_or(&m.species);
        let Some(s) = species.iter().position(|x| x == label) else {
            return Err(RadkinError::Config(format!(
                "seed molecule '{}' (radkin '{label}') is not a dynamic species of the chemistry; add an alias or list it in the geant4.ignore list",
                m.species
            )));
        };
        let pos: Vec3 = std::array::from_fn(|a| m.pos[a] + spec.offset[a]);
        if let Some((lo, hi)) = domain
            && (0..3).any(|a| !(pos[a] > lo[a] && pos[a] < hi[a]))
        {
            seed.clipped += 1;
            continue;
        }
        seed.molecules.push((s as u32, pos));
    }
    if seed.clipped > 0 {
        let msg = format!(
            "{}: {} of {} seed molecules are outside the box",
            spec.file.display(),
            seed.clipped,
            seed.clipped + seed.molecules.len()
        );
        if !policy.clip {
            return Err(RadkinError::Config(format!(
                "{msg}; enlarge the box or set clip = true"
            )));
        }
        log::warn!("{msg}: dropped");
    }
    Ok(seed)
}

/// A run of the IRT stage.
#[derive(Debug, Clone)]
pub struct IrtSimulation {
    pub env: SimEnv,
    pub model: IrtModel,
    pub config: IrtConfig,
    pub seeds: Vec<IrtSeed>,
    /// Start of the stage [s].
    pub t_start: f64,
}

/// Results of an IRT ensemble.
#[derive(Debug, Clone, PartialEq)]
pub struct IrtOutput {
    /// Labels of the dynamic species.
    pub species: Vec<String>,
    pub record_times: Vec<f64>,
    /// Mean particles of each species at each record time.
    pub totals_mean: Vec<Vec<f64>>,
    /// Standard deviation over the trajectories.
    pub totals_std: Vec<Vec<f64>>,
    /// `(time, energy [eV])` of every seed that gives its energy.
    pub seed_energies: Vec<(f64, f64)>,
    pub trajectories: usize,
    pub seed: u64,
    /// Event counts summed over the trajectories.
    pub stats: StageStats,
    /// Reaction strings of the network, by reaction index.
    pub reactions: Vec<String>,
    pub warnings: Vec<String>,
    pub exec_time: f64,
}

impl IrtSimulation {
    /// Load a configuration file with an `[irt]` section.
    pub fn from_toml(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let cfg = Config::load(path)?;
        let irt = IrtConfig::from_toml_file(path)?.ok_or_else(|| {
            RadkinError::Config(format!("{}: missing [irt] section", path.display()))
        })?;
        if !cfg.beam.is_empty() {
            return Err(RadkinError::Config(
                "irt: a configuration with a [beam] cannot run an IRT stage (its molecules come from seeds, not from G-values)".into(),
            ));
        }
        Self::new(SimEnv::from_config(&cfg)?, irt)
    }

    pub fn new(env: SimEnv, config: IrtConfig) -> Result<Self> {
        let mat = SimMatrices::build(&env)?;
        let (table, extra) = config.species_table()?;
        let model = IrtModel::with_diffusion(&env, &mat, &table, &extra)?;
        let domain = match config.boundary {
            BoundaryKind::None => None,
            BoundaryKind::Outflow => config.domain,
        };
        let seeds = config
            .seeds
            .iter()
            .map(|s| load_seed(s, &model.species, &config.policy, domain))
            .collect::<Result<Vec<_>>>()?;
        let first = seeds
            .iter()
            .map(|s| s.spec.time)
            .fold(f64::INFINITY, f64::min);
        let t_start = match config.t_start {
            Some(t) => t,
            None if first.is_finite() => first,
            None => {
                return Err(RadkinError::Config(
                    "irt: no seed and no t_start: nothing to run".into(),
                ));
            }
        };
        if t_start >= config.t_end {
            return Err(RadkinError::Config(format!(
                "irt: t_end ({}) must be above the start time ({t_start})",
                config.t_end
            )));
        }
        if let Some(bad) = seeds
            .iter()
            .find(|s| s.spec.time < t_start || s.spec.time >= config.t_end)
        {
            return Err(RadkinError::Config(format!(
                "irt: the seed {} is injected at {} s, outside the IRT stage [{t_start}, {}) s",
                bad.spec.file.display(),
                bad.spec.time,
                config.t_end
            )));
        }
        Ok(Self {
            env,
            model,
            config,
            seeds,
            t_start,
        })
    }

    /// Record times of the run: the configured ones inside the stage, or
    /// 100 log-spaced times; always with the start and the end.
    pub fn record_times(&self) -> Vec<f64> {
        let (a, b) = (self.t_start, self.config.t_end);
        let mut times = match &self.config.record_times {
            Some(list) => {
                let mut t: Vec<f64> = list.iter().copied().filter(|t| *t > a && *t < b).collect();
                t.push(a);
                t.push(b);
                t
            }
            None => thinned_times(a, b, N_RECORDS),
        };
        times.sort_by(f64::total_cmp);
        times.dedup();
        times
    }

    fn stage_options(&self) -> StageOptions {
        StageOptions {
            t_stop: self.config.t_end,
            start_contact: self.config.start_contact,
            boundary: match (self.config.boundary, self.config.domain) {
                (BoundaryKind::Outflow, Some((lo, hi))) => Boundary::Outflow { lo, hi },
                _ => Boundary::None,
            },
            cutoff_probability: self.config.cutoff_probability,
            molecule_reactions: self.config.background_reactions,
        }
    }

    /// One trajectory: the particle counts at every record time and the
    /// event counts.
    pub fn run_trajectory(&self, seed: u64, stream: u64) -> (Vec<Vec<i64>>, StageStats) {
        let records = self.record_times();
        let mut times = records.clone();
        times.extend(self.seeds.iter().map(|s| s.spec.time));
        times.sort_by(f64::total_cmp);
        times.dedup();
        let mut stage = Stage::new(
            &self.model,
            self.stage_options(),
            SsaRng::new(seed, stream),
            self.t_start,
        );
        let mut totals = Vec::with_capacity(records.len());
        for &t in &times {
            stage.advance_to(t);
            let arriving = self.seeds.iter().filter(|s| s.spec.time == t);
            let mut synced = false;
            for s in arriving {
                if !synced {
                    stage.sync(t);
                    synced = true;
                }
                stage.inject(&s.molecules);
            }
            if t == self.config.t_end {
                stage.sync(t);
            }
            if records.contains(&t) {
                totals.push(stage.counts().to_vec());
            }
        }
        (totals, stage.stats().clone())
    }
}

fn add_stats(sum: &mut StageStats, s: &StageStats) {
    sum.pair_reactions += s.pair_reactions;
    sum.contact_reactions += s.contact_reactions;
    sum.molecule_reactions += s.molecule_reactions;
    sum.passes += s.passes;
    sum.moves += s.moves;
    sum.pairs_sampled += s.pairs_sampled;
    if sum.by_reaction.is_empty() {
        sum.by_reaction = vec![0; s.by_reaction.len()];
        sum.absorbed = vec![0; s.absorbed.len()];
        sum.background_consumed = vec![0; s.background_consumed.len()];
    }
    for (a, b) in sum
        .background_consumed
        .iter_mut()
        .zip(&s.background_consumed)
    {
        *a += b;
    }
    for (a, b) in sum.by_reaction.iter_mut().zip(&s.by_reaction) {
        *a += b;
    }
    for (a, b) in sum.absorbed.iter_mut().zip(&s.absorbed) {
        *a += b;
    }
}

/// Run the ensemble of `sim` with `trajectories` trajectories keyed by
/// `seed`.
pub fn run_irt_with(sim: &IrtSimulation, trajectories: usize, seed: u64) -> Result<IrtOutput> {
    let started = Instant::now();
    if trajectories == 0 {
        return Err(RadkinError::Value(
            "irt: at least one trajectory is needed".into(),
        ));
    }
    let record_times = sim.record_times();
    let n_species = sim.model.n_species();
    let mut moments = vec![vec![Moments::default(); n_species]; record_times.len()];
    let mut stats = StageStats::default();
    let runs = map_ordered(0..trajectories, |k| sim.run_trajectory(seed, k as u64));
    for (totals, s) in &runs {
        for (row, counts) in moments.iter_mut().zip(totals) {
            for (m, c) in row.iter_mut().zip(counts) {
                m.add(*c);
            }
        }
        add_stats(&mut stats, s);
    }
    let of = |f: fn(&Moments) -> f64| -> Vec<Vec<f64>> {
        moments
            .iter()
            .map(|row| row.iter().map(f).collect())
            .collect()
    };
    Ok(IrtOutput {
        species: sim.model.species.clone(),
        totals_mean: of(Moments::mean),
        totals_std: of(Moments::std),
        record_times,
        seed_energies: sim
            .seeds
            .iter()
            .filter_map(|s| s.spec.energy_ev.map(|e| (s.spec.time, e)))
            .collect(),
        trajectories,
        seed,
        stats,
        reactions: sim.env.reaction_strings(),
        warnings: sim.model.warnings.clone(),
        exec_time: started.elapsed().as_secs_f64(),
    })
}

/// Run the ensemble of `sim` with the trajectories and seed of its `[irt]`
/// section (1 and 0 by default).
pub fn run_irt(sim: &IrtSimulation) -> Result<IrtOutput> {
    run_irt_with(
        sim,
        sim.config.trajectories.unwrap_or(1),
        sim.config.seed.unwrap_or(0),
    )
}
