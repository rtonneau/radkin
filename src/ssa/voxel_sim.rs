//! Ensemble of RDME trajectories of a voxel run (`radkin voxel` with
//! `solver = "ssa"`).
//!
//! The run takes the same configuration as the deterministic voxel path:
//! chemistry, `[voxel]` grid, diffusion coefficients and seeds, and
//! `[voxel.run]` times. The splitting and rsoda settings of `[voxel.run]`
//! are not used.
//!
//! Each trajectory starts from a Poisson draw of the initial solution in
//! every voxel, receives the Geant4 seeds as exact integer counts at their
//! times, and is recorded at the snapshot times (full fields) and at the
//! record times (domain totals). Trajectory `i` draws from stream `i` of
//! the seed; the statistics are accumulated in trajectory order.

use std::time::Instant;

use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::ssa::config::{RateCorrection, SsaConfig, hierarchy_depth, suggested_level};
use crate::ssa::ensemble::{Moments, map_ordered};
use crate::ssa::rdme::{Rdme, RdmeSetup, ReversiblePair, is_cubic};
use crate::ssa::rdme_coarsen::dispersion_threshold;
use crate::ssa::rdme_multigrid::{Multigrid, MultigridSetup};
use crate::ssa::rng::SsaRng;
use crate::unit_cell::thinned_times;
use crate::voxel::config::check_divisible;
use crate::voxel::grid::{VoxelGrid, VoxelState};
use crate::voxel::sim::VoxelSimulation;

/// Log-spaced record times of the domain totals, on top of the snapshot and
/// seed times.
pub const N_RECORDS: usize = 100;

/// Trajectories run at once before their sums are accumulated.
const BATCH: usize = 32;

/// Sum and sum of squares of one field value over the trajectories. The
/// sum of squares is a float: it is accumulated in trajectory order, so it
/// is reproducible, and it keeps a snapshot at 16 bytes per value.
#[derive(Debug, Clone, Copy, Default)]
struct FieldMoments {
    sum: i64,
    sum_sq: f64,
}

/// Everything an RDME voxel run produces.
#[derive(Debug, Clone)]
pub struct RdmeOutput {
    /// Dynamic species labels, in state order.
    pub species: Vec<String>,
    /// Times of the snapshots [s].
    pub snapshot_times: Vec<f64>,
    /// Ensemble mean field of every snapshot [µmol/L], on the fine grid.
    pub mean: Vec<VoxelState>,
    /// Standard deviation over the trajectories [µmol/L].
    pub std: Vec<VoxelState>,
    /// Times of the domain totals [s].
    pub record_times: Vec<f64>,
    /// Mean molecules of every species in the domain, per record time.
    pub totals_mean: Vec<Vec<f64>>,
    /// Standard deviation of the totals over the trajectories.
    pub totals_std: Vec<Vec<f64>>,
    /// Mean and standard deviation of the molecules net of the totals at
    /// `t_start` (before any seed), per record time.
    pub net_mean: Vec<Vec<f64>>,
    pub net_std: Vec<Vec<f64>>,
    /// `(time, energy_ev)` of every seed with a deposited energy.
    pub seed_energies: Vec<(f64, f64)>,
    pub trajectories: usize,
    pub seed: u64,
    pub tau_leap: bool,
    /// Mesh correction applied: `"hellander"` or `"none"`.
    pub rate_correction: &'static str,
    /// Events summed over the trajectories.
    pub n_reactions: u64,
    pub n_hops: u64,
    /// `(t, level)` after every change of coarsening level of trajectory 0.
    pub level_history: Vec<(f64, u32)>,
    /// Species level of every species (`[ssa.species_levels]`), state order.
    pub species_levels: Vec<u32>,
    /// Wall time [s].
    pub exec_time: f64,
}

/// What one trajectory returns.
struct TrajectoryOutput {
    /// Fine fields at the snapshot times ([`Rdme::fine_field`] units).
    snapshots: Vec<Vec<i64>>,
    /// Domain totals at the record times.
    totals: Vec<Vec<i64>>,
    /// Totals at `t_start`, before any seed.
    baseline: Vec<i64>,
    n_reactions: u64,
    n_hops: u64,
    level_history: Vec<(f64, u32)>,
}

/// The times a trajectory stops at, in order.
struct Schedule {
    times: Vec<f64>,
    /// Indices into the seed list injected at each time.
    seeds: Vec<Vec<usize>>,
    is_record: Vec<bool>,
    is_snapshot: Vec<bool>,
}

fn schedule(sim: &VoxelSimulation) -> Result<Schedule> {
    let s = &sim.settings;
    if let Some(bad) = sim
        .seeds
        .iter()
        .find(|x| x.spec.time < s.t_start || x.spec.time > s.t_end)
    {
        return Err(RadkinError::Config(format!(
            "seed time {} s is outside the run [{}, {}] s",
            bad.spec.time, s.t_start, s.t_end
        )));
    }
    let snapshots: Vec<f64> = s
        .snapshots
        .iter()
        .copied()
        .filter(|t| *t >= s.t_start && *t <= s.t_end)
        .collect();
    let mut times = thinned_times(s.t_start, s.t_end, N_RECORDS);
    times.extend(&snapshots);
    times.extend(sim.seeds.iter().map(|x| x.spec.time));
    times.sort_by(f64::total_cmp);
    times.dedup();
    let seeds = times
        .iter()
        .map(|t| {
            (0..sim.seeds.len())
                .filter(|i| sim.seeds[*i].spec.time == *t)
                .collect()
        })
        .collect();
    let is_snapshot = times.iter().map(|t| snapshots.contains(t)).collect();
    Ok(Schedule {
        is_record: vec![true; times.len()],
        times,
        seeds,
        is_snapshot,
    })
}

/// `(fine voxel, species, molecules)` entries of every seed: the binned
/// concentrations back to the integer counts they came from.
fn seed_molecules(sim: &VoxelSimulation) -> Result<Vec<Vec<(usize, usize, i64)>>> {
    let unit = sim.voxel.grid.molecule_concentration() * 1e6;
    sim.seeds
        .iter()
        .map(|seed| {
            let st = &seed.amount;
            let mut out = Vec::new();
            for s in 0..st.n_species {
                for (v, amount) in st.field(s).iter().enumerate() {
                    if *amount == 0.0 {
                        continue;
                    }
                    let n = (amount / unit).round();
                    if n < 0.0 || (amount / unit - n).abs() > 1e-6 * n.max(1.0) {
                        return Err(RadkinError::Value(format!(
                            "rdme: seed {} does not hold whole molecules in voxel {v}",
                            seed.spec.file.display()
                        )));
                    }
                    out.push((v, s, n as i64));
                }
            }
            Ok(out)
        })
        .collect()
}

/// Whether the run applies the mesh correction: an unset key means
/// Hellander on a cubic grid and none (with a warning) otherwise; an
/// explicit Hellander on an anisotropic grid is an error.
fn resolve_correction(setting: Option<RateCorrection>, grid: &VoxelGrid) -> Result<bool> {
    let cubic = is_cubic(grid);
    match setting {
        Some(RateCorrection::None) => Ok(false),
        Some(RateCorrection::Hellander) if cubic => Ok(true),
        Some(RateCorrection::Hellander) => Err(RadkinError::Config(format!(
            "ssa.rate_correction = \"hellander\" needs equal voxel spacings, got {:?} nm",
            grid.spacing.map(|h| h * 1e9)
        ))),
        None if cubic => Ok(true),
        None => {
            log::warn!(
                "rdme: the voxel spacings {:?} nm are not equal: the bimolecular propensities \
                 are not mesh-corrected (set ssa.rate_correction = \"none\" to silence this)",
                grid.spacing.map(|h| h * 1e9)
            );
            Ok(false)
        }
    }
}

/// The setup of a run: one level for every species, or species levels.
enum Model {
    Single(RdmeSetup),
    Multi(Box<MultigridSetup>),
}

/// What the ensemble driver needs from one trajectory.
trait VoxelTrajectory {
    fn fill_poisson(&mut self, y: &[f64]);
    fn totals(&self) -> Vec<i64>;
    fn advance_to(&mut self, t: f64);
    fn inject(&mut self, molecules: &[(usize, usize, i64)]);
    fn fine_field(&self) -> Vec<i64>;
    /// Reactions, hops and the coarsening history.
    fn finish(self) -> (u64, u64, Vec<(f64, u32)>);
}

impl VoxelTrajectory for Rdme<'_> {
    fn fill_poisson(&mut self, y: &[f64]) {
        Rdme::fill_poisson(self, y)
    }
    fn totals(&self) -> Vec<i64> {
        Rdme::totals(self)
    }
    fn advance_to(&mut self, t: f64) {
        Rdme::advance_to(self, t)
    }
    fn inject(&mut self, molecules: &[(usize, usize, i64)]) {
        Rdme::inject(self, molecules)
    }
    fn fine_field(&self) -> Vec<i64> {
        Rdme::fine_field(self)
    }
    fn finish(self) -> (u64, u64, Vec<(f64, u32)>) {
        (self.n_reactions, self.n_hops, self.level_history)
    }
}

impl VoxelTrajectory for Multigrid<'_> {
    fn fill_poisson(&mut self, y: &[f64]) {
        Multigrid::fill_poisson(self, y)
    }
    fn totals(&self) -> Vec<i64> {
        Multigrid::totals(self)
    }
    fn advance_to(&mut self, t: f64) {
        Multigrid::advance_to(self, t)
    }
    fn inject(&mut self, molecules: &[(usize, usize, i64)]) {
        Multigrid::inject(self, molecules)
    }
    fn fine_field(&self) -> Vec<i64> {
        Multigrid::fine_field(self)
    }
    fn finish(self) -> (u64, u64, Vec<(f64, u32)>) {
        (self.n_reactions, self.n_hops, self.level_history)
    }
}

/// Run one trajectory through the schedule: Poisson initial solution, then
/// the seeds, records and snapshots at their times.
fn record(
    mut traj: impl VoxelTrajectory,
    plan: &Schedule,
    molecules: &[Vec<(usize, usize, i64)>],
    y0: &[f64],
) -> TrajectoryOutput {
    traj.fill_poisson(y0);
    let baseline = traj.totals();
    let (mut snapshots, mut totals) = (Vec::new(), Vec::new());
    for (k, t) in plan.times.iter().enumerate() {
        traj.advance_to(*t);
        for seed in &plan.seeds[k] {
            traj.inject(&molecules[*seed]);
        }
        if plan.is_record[k] {
            totals.push(traj.totals());
        }
        if plan.is_snapshot[k] {
            snapshots.push(traj.fine_field());
        }
    }
    let (n_reactions, n_hops, level_history) = traj.finish();
    TrajectoryOutput {
        snapshots,
        totals,
        baseline,
        n_reactions,
        n_hops,
        level_history,
    }
}

/// The species level of every dynamic species (`[ssa.species_levels]`),
/// after checking that the grid can hold the hierarchy they need. Logs the
/// suggested level of every background species (non-zero initial
/// concentration, not seeded) from the diffusion coefficients (ADR 0005);
/// nothing is applied.
fn species_levels(
    sim: &VoxelSimulation,
    cfg: &SsaConfig,
    species: &[String],
    coeffs: &[f64],
) -> Result<Vec<u32>> {
    let levels = cfg.resolve_species_levels(species)?;
    let max_coarsen = sim.voxel.max_coarsen_levels;
    let depth = hierarchy_depth(max_coarsen, &levels);
    if depth > max_coarsen {
        let deepest = (0..species.len())
            .find(|s| levels[*s] == depth)
            .unwrap_or(0);
        let key = format!("ssa.species_levels.{}", species[deepest]);
        check_divisible(sim.voxel.grid.dims, depth, &key)?;
    }
    let seeded: Vec<bool> = (0..species.len())
        .map(|s| {
            sim.seeds
                .iter()
                .any(|x| x.amount.field(s).iter().any(|a| *a != 0.0))
        })
        .collect();
    let d_slowest = (0..species.len())
        .filter(|s| seeded[*s] && coeffs[*s] > 0.0)
        .map(|s| coeffs[s])
        .fold(f64::INFINITY, f64::min);
    // The suggestion may use levels the grid does not divide into.
    let max_depth = sim
        .voxel
        .grid
        .dims
        .iter()
        .map(|n| n.trailing_zeros())
        .min()
        .unwrap_or(0);
    let y0 = sim.env.initial_values();
    for (s, name) in species.iter().enumerate() {
        if y0[s] <= 0.0 || seeded[s] {
            continue;
        }
        if let Some(l) = suggested_level(coeffs[s], d_slowest, max_depth) {
            log::info!(
                "rdme: suggested species level of background {name}: {l} \
                 (floor(log2 sqrt(D / D_slowest seeded)), {:.2e} / {:.2e} m^2/s; set {})",
                coeffs[s],
                d_slowest,
                levels[s]
            );
        }
    }
    Ok(levels)
}

/// Log how the mesh correction treats every reversible pair whose
/// bimolecular direction it corrects (ADR 0006), with the factors of the
/// fine grid.
fn log_reversible_pairs(setup: &RdmeSetup, reactions: &[String]) {
    let factor = |j: usize| 1.0 / setup.pole_distance(0, j);
    for p in setup.reversible_pairs() {
        match *p {
            ReversiblePair::Shared {
                bimolecular,
                follower,
            } if setup.pole_distance(0, bimolecular) < 1.0 => log::info!(
                "rdme: reversible pair '{}': both directions take the mesh correction of \
                 '{}' (x{:.3} on the fine grid), so K_eq is kept",
                // The forward direction, as declared.
                reactions[follower.min(bimolecular)],
                reactions[bimolecular],
                factor(bimolecular)
            ),
            ReversiblePair::BothBimolecular { forward, backward }
                if setup.critical_spacing(forward).is_some()
                    || setup.critical_spacing(backward).is_some() =>
            {
                log::warn!(
                    "rdme: reversible pair '{}' is bimolecular in both directions: each keeps \
                     its own mesh correction (x{:.3} forward, x{:.3} backward on the fine \
                     grid), so its K_eq shifts by x{:.3}",
                    reactions[forward],
                    factor(forward),
                    factor(backward),
                    factor(forward) / factor(backward)
                )
            }
            _ => {}
        }
    }
}

/// Run the RDME ensemble of `sim`.
pub fn run_rdme(sim: &VoxelSimulation, cfg: &SsaConfig) -> Result<RdmeOutput> {
    let started = Instant::now();
    let species = sim.species();
    let mat = SimMatrices::build(&sim.env)?;
    let coeffs: Vec<f64> = species.iter().map(|s| sim.voxel.diffusion_of(s)).collect();
    let levels = species_levels(sim, cfg, &species, &coeffs)?;
    let multigrid = levels.iter().any(|l| *l > 0);
    let depth = hierarchy_depth(sim.voxel.max_coarsen_levels, &levels);
    let grid = &sim.voxel.grid;
    let correct = resolve_correction(cfg.rate_correction, grid)?;
    let setup = RdmeSetup::new(&mat, grid, coeffs, depth, correct)?;
    let smallest = grid.spacing.iter().copied().fold(f64::INFINITY, f64::min);
    let reactions = sim.env.reaction_strings();
    if correct {
        log_reversible_pairs(&setup, &reactions);
        for &(j, h) in setup.critical_spacings() {
            if smallest <= h {
                log::warn!(
                    "rdme: reaction '{}' is not mesh-corrected on the fine grid: the voxel \
                     spacing {:.2} nm is at or below its critical spacing {:.2} nm \
                     (0.2527 kappa / D), where no propensity represents it: its rate is \
                     underestimated",
                    reactions[j],
                    smallest * 1e9,
                    h * 1e9,
                );
            }
        }
    } else if let Some((h, j)) = setup.validity_bound()
        && smallest < h
    {
        log::warn!(
            "rdme: the voxel spacing {:.2} nm is below {:.2} nm, the smallest spacing at \
             which reaction '{}' is represented correctly on a grid (0.2527 kappa / D): \
             its rate is underestimated",
            smallest * 1e9,
            h * 1e9,
            reactions[j]
        );
    }
    let z_max = (sim.voxel.max_coarsen_levels > 0)
        .then(|| dispersion_threshold(sim.voxel.coarsen_tol_factor));
    let plan = schedule(sim)?;
    let molecules = seed_molecules(sim)?;
    let y0 = sim.env.initial_values();
    let t_start = sim.settings.t_start;

    let max_levels = setup.max_levels();
    let model = if multigrid {
        let mg = MultigridSetup::new(setup, &levels, sim.voxel.max_coarsen_levels)?;
        Model::Multi(Box::new(mg))
    } else {
        Model::Single(setup)
    };

    let one = |i: usize| -> TrajectoryOutput {
        let rng = SsaRng::new(cfg.seed, i as u64);
        match &model {
            Model::Multi(mg) => {
                let mut traj = Multigrid::new(mg, t_start, rng);
                traj.set_coarsening(z_max);
                traj.set_tau_leaping(cfg.tau_leap.then_some(cfg.tau_epsilon));
                record(traj, &plan, &molecules, &y0)
            }
            Model::Single(setup) => {
                let mut rdme = Rdme::new(setup, t_start, rng);
                rdme.set_coarsening(z_max);
                rdme.set_tau_leaping(cfg.tau_leap.then_some(cfg.tau_epsilon));
                record(rdme, &plan, &molecules, &y0)
            }
        }
    };

    let n_sp = species.len();
    let n_fine = sim.voxel.grid.n_voxels();
    let snapshot_times: Vec<f64> = plan
        .times
        .iter()
        .zip(&plan.is_snapshot)
        .filter(|(_, s)| **s)
        .map(|(t, _)| *t)
        .collect();
    let mut fields = vec![vec![FieldMoments::default(); n_sp * n_fine]; snapshot_times.len()];
    let mut totals = vec![vec![Moments::default(); n_sp]; plan.times.len()];
    let mut net = vec![vec![Moments::default(); n_sp]; plan.times.len()];
    let (mut n_reactions, mut n_hops) = (0u64, 0u64);
    let mut level_history = Vec::new();
    let mut first = 0;
    while first < cfg.trajectories {
        let last = (first + BATCH).min(cfg.trajectories);
        for (i, traj) in map_ordered(first..last, one).into_iter().enumerate() {
            n_reactions += traj.n_reactions;
            n_hops += traj.n_hops;
            if first + i == 0 {
                level_history = traj.level_history;
            }
            for (snap, acc) in traj.snapshots.iter().zip(fields.iter_mut()) {
                for (x, m) in snap.iter().zip(acc.iter_mut()) {
                    m.sum += x;
                    m.sum_sq += (*x as f64) * (*x as f64);
                }
            }
            for (k, row) in traj.totals.iter().enumerate() {
                for (s, x) in row.iter().enumerate() {
                    totals[k][s].add(*x);
                    net[k][s].add(*x - traj.baseline[s]);
                }
            }
        }
        first = last;
    }

    // Field units: 1 / 8^max_levels molecule per fine voxel -> µmol/L.
    let n = cfg.trajectories as f64;
    let unit = sim.voxel.grid.molecule_concentration() * 1e6 / 8f64.powi(max_levels as i32);
    let field = |f: &dyn Fn(&FieldMoments) -> f64| -> Vec<VoxelState> {
        fields
            .iter()
            .map(|acc| VoxelState {
                n_species: n_sp,
                n_voxels: n_fine,
                data: acc.iter().map(|m| f(m) * unit).collect(),
            })
            .collect()
    };
    let mean = field(&|m| m.sum as f64 / n);
    let std = field(&|m| {
        if n < 2.0 {
            return 0.0;
        }
        let s = m.sum as f64;
        ((m.sum_sq - s * s / n).max(0.0) / (n - 1.0)).sqrt()
    });
    let table = |acc: &[Vec<Moments>], f: fn(&Moments) -> f64| -> Vec<Vec<f64>> {
        acc.iter().map(|row| row.iter().map(f).collect()).collect()
    };
    Ok(RdmeOutput {
        species,
        snapshot_times,
        mean,
        std,
        record_times: plan.times,
        totals_mean: table(&totals, Moments::mean),
        totals_std: table(&totals, Moments::std),
        net_mean: table(&net, Moments::mean),
        net_std: table(&net, Moments::std),
        seed_energies: sim
            .seeds
            .iter()
            .filter_map(|s| s.spec.energy_ev.map(|e| (s.spec.time, e)))
            .collect(),
        trajectories: cfg.trajectories,
        seed: cfg.seed,
        tau_leap: cfg.tau_leap,
        rate_correction: if correct { "hellander" } else { "none" },
        n_reactions,
        n_hops,
        level_history,
        species_levels: levels,
        exec_time: started.elapsed().as_secs_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(spacing: [f64; 3]) -> VoxelGrid {
        VoxelGrid::new([2, 2, 2], spacing, [0.0; 3]).unwrap()
    }

    #[test]
    fn correction_defaults_to_hellander_on_cubic_grids_only() {
        let cubic = grid([5e-9; 3]);
        let flat = grid([5e-9, 5e-9, 10e-9]);
        assert!(resolve_correction(None, &cubic).unwrap());
        assert!(!resolve_correction(None, &flat).unwrap());
        assert!(resolve_correction(Some(RateCorrection::Hellander), &cubic).unwrap());
        assert!(!resolve_correction(Some(RateCorrection::None), &cubic).unwrap());
        assert!(!resolve_correction(Some(RateCorrection::None), &flat).unwrap());
        let e = resolve_correction(Some(RateCorrection::Hellander), &flat).unwrap_err();
        assert!(e.to_string().contains("rate_correction"), "{e}");
    }
}
