//! Ensemble of trajectories of a unit cell: mean and standard deviation of
//! the molecule counts on the output time grid.
//!
//! Trajectory `i` draws from stream `i` of the seed, whatever thread runs
//! it, and the statistics are exact integer sums accumulated in trajectory
//! order: the result does not depend on the number of threads.

use std::time::Instant;

use crate::chemistry::env::SimEnv;
use crate::error::{RadkinError, Result};
use crate::physics::Beam;
use crate::ssa::config::SsaConfig;
use crate::ssa::exact::{Exact, Stepper, Trajectory, initial_counts, run_trajectory};
use crate::ssa::model::SsaModel;
use crate::ssa::rng::SsaRng;
use crate::ssa::schedule::beam_segments;
use crate::unit_cell::{UnitCell, thinned_times};

/// Number of log-spaced output times when none is given.
pub const DEFAULT_N_POINTS: usize = 200;

/// Trajectories run at once before their sums are accumulated.
const BATCH: usize = 256;

/// `f(i)` for `i` in `range`, results in index order.
#[cfg(feature = "parallel")]
pub(crate) fn map_ordered<T, F>(range: std::ops::Range<usize>, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    use rayon::prelude::*;
    range.into_par_iter().map(f).collect()
}

/// `f(i)` for `i` in `range`, results in index order.
#[cfg(not(feature = "parallel"))]
pub(crate) fn map_ordered<T, F>(range: std::ops::Range<usize>, f: F) -> Vec<T>
where
    F: Fn(usize) -> T,
{
    range.map(f).collect()
}

/// Exact sums of integer samples: mean and sample standard deviation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct Moments {
    n: u64,
    sum: i128,
    sum_sq: u128,
}

impl Moments {
    pub(crate) fn add(&mut self, x: i64) {
        self.n += 1;
        self.sum += i128::from(x);
        self.sum_sq += (i128::from(x) * i128::from(x)) as u128;
    }

    pub(crate) fn mean(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        self.sum as f64 / self.n as f64
    }

    /// Sample standard deviation (0 for fewer than two samples).
    pub(crate) fn std(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        // n * sum_sq - sum^2 >= 0 exactly (Cauchy-Schwarz on integers).
        let n = u128::from(self.n);
        let num = n * self.sum_sq - (self.sum * self.sum) as u128;
        (num as f64 / (self.n as f64 * (self.n - 1) as f64)).sqrt()
    }
}

/// Result of a unit-cell ensemble.
#[derive(Debug, Clone)]
pub struct SsaResult {
    pub env: SimEnv,
    pub beam: Option<Beam>,
    pub let_: Option<f64>,
    /// Output times [s].
    pub t: Vec<f64>,
    /// Ensemble mean of the dynamic species [µmol/L], `mean[i][s]` at `t[i]`.
    pub mean: Vec<Vec<f64>>,
    /// Sample standard deviation over the trajectories [µmol/L].
    pub std: Vec<Vec<f64>>,
    /// Volume [L].
    pub volume: f64,
    pub trajectories: usize,
    pub seed: u64,
    pub tau_leap: bool,
    /// Events fired, summed over the trajectories.
    pub n_events: u64,
    /// Wall time [s].
    pub exec_time: f64,
}

impl SsaResult {
    /// Standard error of the mean [µmol/L].
    pub fn sem(&self) -> Vec<Vec<f64>> {
        let n = (self.trajectories as f64).sqrt();
        self.std
            .iter()
            .map(|row| row.iter().map(|s| s / n).collect())
            .collect()
    }

    /// One molecule in the volume [µmol/L].
    pub fn molecule(&self) -> f64 {
        1e6 / (crate::consts::AVOGADRO * self.volume)
    }

    /// Species labels (dynamic then fixed) and rows of mean concentrations
    /// [µmol/L], as [`crate::result::OdeResult::concentrations`].
    pub fn concentrations(&self) -> (Vec<String>, Vec<Vec<f64>>) {
        let mut labels = self.env.species.sorted_dyn_labels();
        let cst: Vec<f64> = self
            .env
            .species
            .sorted_cst()
            .iter()
            .map(|s| s.initial_cc(self.env.temperature) * 1e6)
            .collect();
        labels.extend(self.env.species.sorted_cst_labels());
        let rows = self
            .mean
            .iter()
            .map(|y| {
                let mut r = y.clone();
                r.extend(&cst);
                r
            })
            .collect();
        (labels, rows)
    }
}

/// What an ensemble needs beyond the `[ssa]` table.
#[derive(Debug, Clone, Default)]
pub struct EnsembleOptions {
    /// Output times; `None`: [`DEFAULT_N_POINTS`] log-spaced times.
    pub t_eval: Option<Vec<f64>>,
    /// Initial state [µmol/L]; `None`: the configuration's.
    pub y0: Option<Vec<f64>>,
    /// Ignore the beam (no source event).
    pub no_beam: bool,
}

/// Stepper of one trajectory for the `[ssa]` settings.
pub(crate) fn make_stepper(cfg: &SsaConfig) -> Result<Box<dyn Stepper>> {
    if cfg.tau_leap {
        return Err(RadkinError::Config(
            "ssa.tau_leap = true is not available".into(),
        ));
    }
    Ok(Box::new(Exact))
}

/// Run the ensemble of `cell` over `t_span`.
///
/// The current beam of the cell (if any) gives the sources; a cell without
/// beam runs the reactions only.
pub fn run_ensemble(
    cell: &UnitCell,
    cfg: &SsaConfig,
    t_span: (f64, f64),
    options: &EnsembleOptions,
) -> Result<SsaResult> {
    let started = Instant::now();
    let volume = cfg.require_volume()?;
    let model = SsaModel::new(&cell.mat, volume)?;
    let beam = if options.no_beam || cell.beams.is_empty() {
        None
    } else {
        Some(cell.beams.current()?.clone())
    };
    let segments = beam_segments(beam.as_ref(), t_span)?;
    let t_eval = match &options.t_eval {
        Some(t) => t.clone(),
        None => thinned_times(t_span.0, t_span.1, DEFAULT_N_POINTS),
    };
    crate::ssa::exact::check_schedule(&segments, &t_eval)?;
    let y0 = match &options.y0 {
        Some(y) if y.len() == model.n_dyn => y.clone(),
        Some(y) => {
            return Err(RadkinError::Value(format!(
                "ssa: y0 has {} values for {} dynamic species",
                y.len(),
                model.n_dyn
            )));
        }
        None => cell.env.initial_values(),
    };
    make_stepper(cfg)?;

    let one = |i: usize| -> Result<Trajectory> {
        let mut rng = SsaRng::new(cfg.seed, i as u64);
        let counts = initial_counts(&model, &y0, &mut rng);
        let mut stepper = make_stepper(cfg)?;
        run_trajectory(
            &model,
            &segments,
            counts,
            &t_eval,
            stepper.as_mut(),
            &mut rng,
        )
    };
    let mut moments = vec![vec![Moments::default(); model.n_dyn]; t_eval.len()];
    let mut n_events = 0u64;
    let mut first = 0;
    while first < cfg.trajectories {
        let last = (first + BATCH).min(cfg.trajectories);
        for traj in map_ordered(first..last, one) {
            let traj = traj?;
            n_events += traj.n_events;
            for (row, m) in traj.rows.iter().zip(moments.iter_mut()) {
                for (x, acc) in row.iter().zip(m.iter_mut()) {
                    acc.add(*x);
                }
            }
        }
        first = last;
    }
    let k = model.molecules_per_umol();
    let table = |f: fn(&Moments) -> f64| -> Vec<Vec<f64>> {
        moments
            .iter()
            .map(|row| row.iter().map(|m| f(m) / k).collect())
            .collect()
    };
    Ok(SsaResult {
        env: cell.env.clone(),
        let_: beam.as_ref().and(cell.beams.current_let_opt()),
        beam,
        t: t_eval,
        mean: table(Moments::mean),
        std: table(Moments::std),
        volume,
        trajectories: cfg.trajectories,
        seed: cfg.seed,
        tau_leap: cfg.tau_leap,
        n_events,
        exec_time: started.elapsed().as_secs_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn moments_are_exact() {
        let mut m = Moments::default();
        assert_eq!((m.mean(), m.std()), (0.0, 0.0));
        for x in [2, 4, 4, 4, 5, 5, 7, 9] {
            m.add(x);
        }
        assert_eq!(m.mean(), 5.0);
        // Sample variance 32 / 7.
        assert!((m.std() - (32.0f64 / 7.0).sqrt()).abs() < 1e-15);
        let mut one = Moments::default();
        one.add(3);
        assert_eq!((one.mean(), one.std()), (3.0, 0.0));
        let mut big = Moments::default();
        for _ in 0..1000 {
            big.add(4_000_000_000);
        }
        assert_eq!((big.mean(), big.std()), (4e9, 0.0));
    }

    #[test]
    fn map_ordered_keeps_the_index_order() {
        let v = map_ordered(0..1000, |i| i * i);
        assert!(v.iter().enumerate().all(|(i, x)| *x == i * i));
    }
}
