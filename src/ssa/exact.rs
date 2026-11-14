//! Exact stochastic simulation of one well-mixed volume (Gillespie's direct
//! method).
//!
//! The dose rate is constant over a [`Segment`]. When a segment ends, the
//! pending event is discarded and a new one is drawn with the propensities
//! of the next segment: waiting times are memoryless, so this is exact for
//! piecewise-constant propensities.

use crate::error::{RadkinError, Result};
use crate::ssa::model::SsaModel;
use crate::ssa::rng::SsaRng;

/// A time span at a constant dose rate [Gy/s].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub start: f64,
    pub stop: f64,
    pub dose_rate: f64,
}

/// State of one trajectory.
#[derive(Debug, Clone)]
pub struct Cell<'a> {
    pub model: &'a SsaModel,
    /// Molecules of every dynamic species.
    pub counts: Vec<i64>,
    pub t: f64,
    /// Reaction and source events fired so far.
    pub n_events: u64,
    /// Propensity of every reaction for `counts`.
    pub(crate) a: Vec<f64>,
}

impl<'a> Cell<'a> {
    pub fn new(model: &'a SsaModel, counts: Vec<i64>, t: f64) -> Self {
        assert_eq!(counts.len(), model.n_dyn, "ssa: counts length");
        let mut cell = Self {
            model,
            a: vec![0.0; model.reactions.len()],
            counts,
            t,
            n_events: 0,
        };
        cell.refresh();
        cell
    }

    /// Recompute every propensity from the counts.
    pub(crate) fn refresh(&mut self) {
        for j in 0..self.a.len() {
            self.a[j] = self.model.propensity(j, &self.counts);
        }
    }

    pub(crate) fn fire_reaction(&mut self, j: usize) {
        let model = self.model;
        for &(s, d) in &model.reactions[j].delta {
            self.counts[s] += d;
            debug_assert!(self.counts[s] >= 0, "ssa: negative count");
        }
        for &r in &model.affects[j] {
            self.a[r] = model.propensity(r, &self.counts);
        }
        self.n_events += 1;
    }

    pub(crate) fn fire_source(&mut self, i: usize) {
        let model = self.model;
        self.counts[model.sources[i].species] += 1;
        for &r in &model.source_affects[i] {
            self.a[r] = model.propensity(r, &self.counts);
        }
        self.n_events += 1;
    }

    /// One exact event before `t_stop` at the dose rate `dose_rate`, or the
    /// time set to `t_stop` if none happens before it. Returns whether an
    /// event fired.
    pub(crate) fn step(&mut self, t_stop: f64, dose_rate: f64, rng: &mut SsaRng) -> bool {
        let model = self.model;
        let a_reac: f64 = self.a.iter().sum();
        let a_src: f64 = if dose_rate > 0.0 {
            model.sources.iter().map(|s| s.per_gray * dose_rate).sum()
        } else {
            0.0
        };
        let a0 = a_reac + a_src;
        if a0 <= 0.0 {
            self.t = t_stop;
            return false;
        }
        let t_next = self.t + rng.exp1() / a0;
        if t_next >= t_stop {
            self.t = t_stop;
            return false;
        }
        self.t = t_next;
        let target = rng.uniform() * a0;
        let mut acc = 0.0;
        if a_src > 0.0 {
            for (i, s) in model.sources.iter().enumerate() {
                acc += s.per_gray * dose_rate;
                if target < acc {
                    self.fire_source(i);
                    return true;
                }
            }
        }
        let mut last = None;
        for j in 0..self.a.len() {
            let a = self.a[j];
            if a > 0.0 {
                acc += a;
                last = Some(j);
                if target < acc {
                    self.fire_reaction(j);
                    return true;
                }
            }
        }
        // Round-off left the target past the last partial sum.
        match last {
            Some(j) => self.fire_reaction(j),
            None => self.fire_source(model.sources.len() - 1),
        }
        true
    }

    /// Advance to `t_stop` at the constant dose rate `dose_rate` [Gy/s].
    pub fn advance(&mut self, t_stop: f64, dose_rate: f64, rng: &mut SsaRng) {
        while self.t < t_stop {
            self.step(t_stop, dose_rate, rng);
        }
    }
}

/// Counts of one trajectory at the output times.
#[derive(Debug, Clone, PartialEq)]
pub struct Trajectory {
    /// `rows[i][s]`: molecules of species `s` at `t_eval[i]`.
    pub rows: Vec<Vec<i64>>,
    pub n_events: u64,
}

/// Check that the segments are contiguous and that the output times are
/// sorted and inside them.
pub fn check_schedule(segments: &[Segment], t_eval: &[f64]) -> Result<()> {
    let (Some(first), Some(last)) = (segments.first(), segments.last()) else {
        return Err(RadkinError::Value("ssa: no time segment".into()));
    };
    for s in segments {
        if s.start.is_nan() || s.stop.is_nan() || s.stop < s.start || s.dose_rate.is_nan() {
            return Err(RadkinError::Value(format!("ssa: invalid segment {s:?}")));
        }
    }
    if segments.windows(2).any(|w| w[0].stop != w[1].start) {
        return Err(RadkinError::Value(
            "ssa: time segments must be contiguous".into(),
        ));
    }
    if t_eval.windows(2).any(|w| w[1].is_nan() || w[1] < w[0]) {
        return Err(RadkinError::Value(
            "ssa: output times must be sorted".into(),
        ));
    }
    if let (Some(a), Some(b)) = (t_eval.first(), t_eval.last())
        && (a.is_nan() || *a < first.start || *b > last.stop)
    {
        return Err(RadkinError::Value(format!(
            "ssa: output times must lie in [{}, {}] s",
            first.start, last.stop
        )));
    }
    Ok(())
}

/// How a trajectory is stepped.
pub trait Stepper {
    /// Advance `cell` to `t_stop` at a constant dose rate.
    fn advance(&mut self, cell: &mut Cell<'_>, t_stop: f64, dose_rate: f64, rng: &mut SsaRng);
}

/// Exact direct method.
#[derive(Debug, Clone, Copy, Default)]
pub struct Exact;

impl Stepper for Exact {
    fn advance(&mut self, cell: &mut Cell<'_>, t_stop: f64, dose_rate: f64, rng: &mut SsaRng) {
        cell.advance(t_stop, dose_rate, rng);
    }
}

/// One trajectory from `counts0` at `segments[0].start`, recorded at `t_eval`.
pub fn run_trajectory(
    model: &SsaModel,
    segments: &[Segment],
    counts0: Vec<i64>,
    t_eval: &[f64],
    stepper: &mut dyn Stepper,
    rng: &mut SsaRng,
) -> Result<Trajectory> {
    check_schedule(segments, t_eval)?;
    let mut cell = Cell::new(model, counts0, segments[0].start);
    let mut rows = Vec::with_capacity(t_eval.len());
    let mut seg = 0;
    for &t_out in t_eval {
        while cell.t < t_out {
            let s = &segments[seg];
            stepper.advance(&mut cell, s.stop.min(t_out), s.dose_rate, rng);
            if cell.t >= s.stop && seg + 1 < segments.len() {
                seg += 1;
            }
        }
        rows.push(cell.counts.clone());
    }
    Ok(Trajectory {
        rows,
        n_events: cell.n_events,
    })
}

/// Molecule counts of the state `y` [µmol/L]: the integer part of
/// `c * N_A * V`, plus one with the probability of its fractional part, so
/// the mean is exact and no more variance than needed is added.
pub fn initial_counts(model: &SsaModel, y: &[f64], rng: &mut SsaRng) -> Vec<i64> {
    let k = model.molecules_per_umol();
    y.iter()
        .map(|v| {
            let x = (v * k).max(0.0);
            let whole = x.floor();
            whole as i64 + i64::from(rng.uniform() < x - whole)
        })
        .collect()
}
