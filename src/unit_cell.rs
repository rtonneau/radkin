//! Simulation driver (port of radiopyo `simulation/unit_cell.py`).

use std::path::Path;
use std::time::Instant;

use crate::chemistry::env::SimEnv;
use crate::config::Config;
use crate::consts::{
    DEFAULT_ATOL, DEFAULT_RTOL, DEFAULT_TRANS_TIME, FIRST_STEP, MIN_BEAM_ON_TIME, T_ZERO,
    THINNING_FIRST_OFFSET, data_dir,
};
use crate::error::{RadkinError, Result};
use crate::gvalues::GValuesBoscolo;
use crate::matrices::SimMatrices;
use crate::physics::{Beam, BeamCollection};
use crate::result::{OdeResult, ShortOnSpan};
use crate::solver::{DoseRateSource, Solution, SolveOptions, solve};

/// Solver options of a run (radiopyo `RunOptions`).
#[derive(Debug, Clone, PartialEq)]
pub struct RunOptions {
    pub atol: Option<f64>,
    pub rtol: Option<f64>,
    pub max_step: Option<f64>,
    pub first_step: Option<f64>,
    pub use_beam: bool,
    pub use_jac: bool,
    pub t_eval: Option<Vec<f64>>,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            atol: None,
            rtol: None,
            max_step: None,
            first_step: None,
            use_beam: true,
            use_jac: true,
            t_eval: None,
        }
    }
}

/// Options of chunked and single-pulse runs.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ChunkOptions {
    /// `atol`, `rtol`, `first_step`, `use_beam`, `use_jac` (max_step and t_eval are per chunk).
    pub run: RunOptions,
    pub max_step_size_on: Option<f64>,
    pub max_step_size_off: Option<f64>,
    pub max_step_size_trans: Option<f64>,
    pub max_step_size_final: Option<f64>,
    pub trans_time: Option<f64>,
    pub points_per_chunk_on: Option<usize>,
    pub points_per_chunk_off: Option<usize>,
    pub y0: Option<Vec<f64>>,
}

/// One chunk of a chunked run, at a constant dose rate.
#[derive(Debug, Clone, PartialEq)]
pub struct SimChunk {
    pub start: f64,
    pub stop: f64,
    pub max_step_size: Option<f64>,
    pub beam_on: bool,
    pub dose_rate: Option<f64>,
    pub n_points: Option<usize>,
}

/// numpy `linspace(start, stop, num)` (endpoint included).
pub fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    if num == 0 {
        return Vec::new();
    }
    if num == 1 {
        return vec![start];
    }
    let div = (num - 1) as f64;
    let delta = stop - start;
    let step = delta / div;
    let mut y: Vec<f64> = (0..num)
        .map(|i| {
            if step == 0.0 {
                i as f64 / div * delta + start
            } else {
                i as f64 * step + start
            }
        })
        .collect();
    y[num - 1] = stop;
    y
}

/// numpy `geomspace(start, stop, num)`.
pub fn geomspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    let mut out: Vec<f64> = linspace(start.log10(), stop.log10(), num)
        .into_iter()
        .map(|e| 10f64.powf(e))
        .collect();
    if num > 0 {
        out[0] = start;
    }
    if num > 1 {
        out[num - 1] = stop;
    }
    out
}

fn unique_sorted(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(f64::total_cmp);
    v.dedup();
    v
}

/// `n_points` log-spaced times over `[start, stop]`, both kept (radiopyo `thinned_times`).
pub fn thinned_times(start: f64, stop: f64, n_points: usize) -> Vec<f64> {
    let span = stop - start;
    if span.is_nan() || span <= 0.0 {
        return unique_sorted(vec![start, stop]);
    }
    let offsets = geomspace(span * THINNING_FIRST_OFFSET, span, n_points - 1);
    let mut times = vec![start];
    times.extend(offsets.iter().map(|o| start + o));
    let last = times.len() - 1;
    times[last] = stop;
    unique_sorted(times.into_iter().map(|t| t.clamp(start, stop)).collect())
}

/// An empty span `(t, t)` becomes `(t/100, t)`.
pub fn widen_empty_span(t_span: (f64, f64)) -> (f64, f64) {
    if t_span.0 == t_span.1 {
        (t_span.0 / 100.0, t_span.1)
    } else {
        t_span
    }
}

/// A run-time ON span shorter than `min_on_time` (with roundoff slack).
pub fn is_short_on_span(t_span: (f64, f64), min_on_time: f64) -> bool {
    let slack = min_on_time * 1e-9 + 4.0 * f64::EPSILON * t_span.0.abs().max(t_span.1.abs());
    t_span.1 - t_span.0 < min_on_time - slack
}

/// Dose rate of a chunk: the beam at its midpoint.
pub fn chunk_dose_rate(beam: &Beam, start: f64, stop: f64) -> f64 {
    let t = if stop.is_finite() {
        0.5 * (start + stop)
    } else {
        start
    };
    beam.at(t)
}

/// `None` or an integer >= 2.
pub fn check_points_per_chunk(value: Option<usize>, name: &str) -> Result<()> {
    match value {
        Some(v) if v < 2 => Err(RadkinError::Value(format!(
            "{name} must be None or an integer >= 2, got {v}"
        ))),
        _ => Ok(()),
    }
}

fn make_chunks(
    boundaries: &[f64],
    beam_on: &[bool],
    max_steps: &[Option<f64>],
    beam: &Beam,
    n_points: &[Option<usize>],
) -> Vec<SimChunk> {
    (0..beam_on.len())
        .map(|i| SimChunk {
            start: boundaries[i],
            stop: boundaries[i + 1],
            max_step_size: max_steps[i],
            beam_on: beam_on[i],
            dose_rate: Some(chunk_dose_rate(beam, boundaries[i], boundaries[i + 1])),
            n_points: n_points[i],
        })
        .collect()
}

/// A homogeneous cell: chemistry, matrices, beams and G-value data.
#[derive(Debug, Clone)]
pub struct UnitCell {
    pub env: SimEnv,
    pub mat: SimMatrices,
    pub beams: BeamCollection,
    pub gvalues: GValuesBoscolo,
    pub last_burst_dose_rate: Option<f64>,
    pub last_short_on_span: Option<ShortOnSpan>,
    pub last_exec_time: f64,
}

impl UnitCell {
    /// Cell from a resolved configuration; the config's current beam is set
    /// with its LET (an error if it has none, as radiopyo).
    pub fn from_config(cfg: &Config) -> Result<Self> {
        let env = SimEnv::from_config(cfg)?;
        let mat = SimMatrices::build(&env)?;
        let gvalues = GValuesBoscolo::load(&data_dir())?;
        let beams = BeamCollection::from_config(&cfg.beam)?;
        let mut cell = Self {
            env,
            mat,
            beams: BeamCollection::new(),
            gvalues,
            last_burst_dose_rate: None,
            last_short_on_span: None,
            last_exec_time: -1.0,
        };
        if !beams.is_empty() {
            let name = beams.current_name().to_string();
            cell.beams = beams;
            cell.use_beam(&name, None)?;
        }
        Ok(cell)
    }

    /// Load a TOML file (includes resolved).
    pub fn from_toml(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_config(&Config::load(path)?)
    }

    fn update_matrices(&mut self) -> Result<()> {
        self.mat = SimMatrices::build(&self.env)?;
        Ok(())
    }

    /// Current beam.
    pub fn beam(&self) -> Result<&Beam> {
        self.beams.current()
    }

    /// Set the temperature [degC] and rebuild the matrices.
    pub fn set_temperature(&mut self, celsius: f64) -> Result<()> {
        self.env.temperature = celsius;
        self.update_matrices()
    }

    /// Add `beam` under `label`, make it current and apply `let_` (or the
    /// LET of the previous current beam).
    pub fn set_beam(&mut self, beam: Beam, label: &str, let_: Option<f64>) -> Result<()> {
        let let_ = match let_ {
            Some(l) => l,
            None => self.beams.current_let()?,
        };
        self.beams.add_beam(beam, label, Some(let_), true);
        self.change_let(let_)
    }

    /// Make the beam `label` current; `let_` overrides its LET.
    pub fn use_beam(&mut self, label: &str, let_: Option<f64>) -> Result<()> {
        let let_ = match let_ {
            Some(l) => Some(l),
            None => self.beams.current_let_opt(),
        };
        self.beams.use_beam(label)?;
        let let_ = match let_ {
            Some(l) => l,
            None => self.beams.current_let()?,
        };
        self.change_let(let_)
    }

    /// Rewrite the LET-dependent G-values and rebuild the matrices.
    pub fn change_let(&mut self, let_: f64) -> Result<()> {
        self.beams.update_let(Some(let_), None)?;
        let peak = self.beams.current()?.dose_rate_peak;
        self.env.apply_let(&self.gvalues, let_, Some(peak))?;
        self.update_matrices()
    }

    fn require_beam(&self) -> Result<()> {
        if self.beams.is_empty() {
            Err(RadkinError::NoBeamDefined)
        } else {
            Ok(())
        }
    }

    fn add_radiolytic_yield(&self, y: &mut [f64], duration: f64, dose_rate: f64) -> Result<()> {
        for r in &self.env.reactions.radiolytic {
            if self.env.species.is_cst(&r.species) {
                continue;
            }
            y[self.env.species.index_of(&r.species)?] +=
                r.initial_cc_before(duration, dose_rate) * 1e6;
        }
        Ok(())
    }

    fn solve_options(&self, t_span: (f64, f64), o: &RunOptions) -> SolveOptions {
        let mut first_step = o.first_step.unwrap_or(FIRST_STEP);
        if let Some(m) = o.max_step
            && first_step >= m
        {
            first_step = m / 10.0;
        }
        if first_step > t_span.1 - t_span.0 {
            first_step = (t_span.1 - t_span.0) / 1e3;
        }
        SolveOptions {
            rtol: o.rtol.unwrap_or(DEFAULT_RTOL),
            atol: o.atol.unwrap_or(DEFAULT_ATOL),
            first_step,
            max_step: o.max_step,
            use_jac: o.use_jac,
            t_eval: o.t_eval.clone(),
        }
    }

    /// One span (radiopyo `_raw_run`).
    fn raw_run(
        &mut self,
        t_span: (f64, f64),
        options: &RunOptions,
        y0: Option<&[f64]>,
        beam_on: bool,
        dose_rate: Option<f64>,
        n_points: Option<usize>,
    ) -> Result<Solution> {
        self.require_beam()?;
        let mut t_span = widen_empty_span(t_span);
        let mut solve_opts = self.solve_options(t_span, options);
        let mut y_start = match y0 {
            Some(y) => y.to_vec(),
            None => self.env.initial_values(),
        };
        let beam = self.beams.current()?.clone();
        let source = match (options.use_beam, dose_rate) {
            (false, _) => DoseRateSource::None,
            (true, None) => DoseRateSource::Beam(beam.clone()),
            (true, Some(d)) => DoseRateSource::Constant(d),
        };
        self.last_burst_dose_rate = None;
        self.last_short_on_span = None;
        if beam_on && options.use_beam {
            if is_short_on_span(t_span, MIN_BEAM_ON_TIME) {
                return self.short_on_span_step(t_span, y_start, dose_rate, options, &beam);
            }
            t_span = (t_span.0 + T_ZERO, t_span.1);
            if t_span.0 >= t_span.1 {
                return Err(RadkinError::Value(format!(
                    "Too Short beam pulse ({}, {})",
                    t_span.0, t_span.1
                )));
            }
            let dr = dose_rate.unwrap_or_else(|| beam.at(t_span.0));
            self.last_burst_dose_rate = Some(dr);
            self.add_radiolytic_yield(&mut y_start, T_ZERO, dr)?;
            if let Some(te) = &solve_opts.t_eval {
                let kept: Vec<f64> = te.iter().copied().filter(|&t| t >= t_span.0).collect();
                if kept.len() != te.len() {
                    log::warn!(
                        "{} t_eval time(s) before {} s (end of the skipped initial burst) dropped",
                        te.len() - kept.len(),
                        t_span.0
                    );
                    solve_opts.t_eval = Some(kept);
                }
            }
        }
        if let Some(n) = n_points
            && t_span.1.is_finite()
        {
            solve_opts.t_eval = Some(thinned_times(t_span.0, t_span.1, n));
        }
        let t0 = Instant::now();
        let sol = solve(&self.mat, &source, t_span, &y_start, &solve_opts)?;
        self.last_exec_time = t0.elapsed().as_secs_f64();
        Ok(sol)
    }

    fn short_on_span_step(
        &mut self,
        t_span: (f64, f64),
        mut y: Vec<f64>,
        dose_rate: Option<f64>,
        options: &RunOptions,
        beam: &Beam,
    ) -> Result<Solution> {
        let (start, stop) = t_span;
        let dr = dose_rate.unwrap_or_else(|| chunk_dose_rate(beam, start, stop));
        let duration = stop - start;
        self.add_radiolytic_yield(&mut y, duration, dr)?;
        self.last_short_on_span = Some(ShortOnSpan {
            start,
            stop,
            dose_rate: dr,
        });
        self.last_exec_time = 0.0;
        log::warn!(
            "Beam ON span [{start:.6e}, {stop:.6e}] s ({duration:.3e} s) is shorter than MIN_BEAM_ON_TIME={MIN_BEAM_ON_TIME:.3e} s: not integrated, its radiolytic yield at {dr:.3e} Gy/s is added analytically"
        );
        if options.t_eval.as_ref().is_some_and(|t| !t.is_empty()) {
            log::warn!("t_eval ignored on the short ON span: one point stored at {stop:.6e} s");
        }
        Ok(Solution {
            t: vec![stop],
            y: vec![y],
            ..Default::default()
        })
    }

    fn new_result(&self) -> Result<OdeResult> {
        Ok(OdeResult::new(
            self.env.clone(),
            self.mat.clone(),
            self.beams.current()?.clone(),
            self.beams.current_let_opt(),
            T_ZERO,
        ))
    }

    fn last_bursts(&self) -> Vec<f64> {
        self.last_burst_dose_rate.into_iter().collect()
    }

    /// Single run over `t_span`, the beam read at every evaluation.
    pub fn run(
        &mut self,
        t_span: (f64, f64),
        options: &RunOptions,
        y0: Option<&[f64]>,
    ) -> Result<OdeResult> {
        self.require_beam()?;
        let sol = self.raw_run(t_span, options, y0, true, None, None)?;
        let mut res = self.new_result()?;
        res.append(
            sol,
            self.last_exec_time,
            &self.last_bursts(),
            None,
            self.last_short_on_span,
            false,
        );
        Ok(res)
    }

    /// Run from `time_start` to the time of the maximum dose.
    pub fn run_till_max_dose(
        &mut self,
        time_start: f64,
        options: &RunOptions,
    ) -> Result<OdeResult> {
        let t_max = self.beam()?.time_at_max_dose;
        if t_max.is_infinite() {
            return Err(RadkinError::Value(
                "Unable to use 'run_till_max_dose', time at max dose is inf".into(),
            ));
        }
        self.run((time_start, t_max), options, None)
    }

    /// Chunks cut at the beam transitions (radiopyo `prepare_chunked_run`).
    pub fn prepare_chunked_run(
        &self,
        t_span: (f64, f64),
        o: &ChunkOptions,
    ) -> Result<Vec<SimChunk>> {
        self.require_beam()?;
        check_points_per_chunk(o.points_per_chunk_on, "points_per_chunk_on")?;
        check_points_per_chunk(o.points_per_chunk_off, "points_per_chunk_off")?;
        let beam = self.beam()?;
        let intervals = beam.make_time_intervals(t_span.1, t_span.0);
        let n = intervals.len();
        let beam_on: Vec<bool> = intervals[..n - 1].iter().map(|r| r[1] > 0.0).collect();
        let mut max_steps: Vec<Option<f64>> = beam_on[..beam_on.len() - 1]
            .iter()
            .map(|&on| {
                if on {
                    o.max_step_size_on
                } else {
                    o.max_step_size_off
                }
            })
            .collect();
        let last = o.max_step_size_final.or(if intervals[n - 2][1] > 0.0 {
            o.max_step_size_on
        } else {
            o.max_step_size_off
        });
        max_steps.push(last);
        let boundaries: Vec<f64> = intervals.iter().map(|r| r[0]).collect();
        let points: Vec<Option<usize>> = beam_on
            .iter()
            .map(|&on| {
                if on {
                    o.points_per_chunk_on
                } else {
                    o.points_per_chunk_off
                }
            })
            .collect();
        Ok(make_chunks(
            &boundaries,
            &beam_on,
            &max_steps,
            beam,
            &points,
        ))
    }

    /// Chunks of a single-pulse run (radiopyo `prepare_single_pulse`).
    pub fn prepare_single_pulse(
        &self,
        t_span: (f64, f64),
        o: &ChunkOptions,
    ) -> Result<Vec<SimChunk>> {
        self.require_beam()?;
        check_points_per_chunk(o.points_per_chunk_on, "points_per_chunk_on")?;
        check_points_per_chunk(o.points_per_chunk_off, "points_per_chunk_off")?;
        let beam = self.beam()?;
        let trans_time = o.trans_time.unwrap_or(DEFAULT_TRANS_TIME);
        let t_max = beam.time_at_max_dose;
        let t_beam_off = t_max + trans_time;
        let (t_start, t_stop) = t_span;
        let mut boundaries = vec![t_start];
        let mut beam_on: Vec<bool> = Vec::new();
        let mut max_steps: Vec<Option<f64>> = Vec::new();
        let t_pre = t_max.min(t_stop);
        if t_pre > t_start {
            let intervals = beam.make_time_intervals(t_pre, t_start);
            boundaries = intervals.iter().map(|r| r[0]).collect();
            beam_on = intervals[..intervals.len() - 1]
                .iter()
                .map(|r| r[1] > 0.0)
                .collect();
            max_steps = beam_on
                .iter()
                .map(|&on| {
                    if on {
                        o.max_step_size_on
                    } else {
                        o.max_step_size_trans
                    }
                })
                .collect();
        }
        for (stop, max_step) in [
            (t_beam_off, o.max_step_size_trans),
            (t_stop, o.max_step_size_final),
        ] {
            let stop = stop.min(t_stop);
            if stop > *boundaries.last().unwrap() {
                boundaries.push(stop);
                beam_on.push(false);
                max_steps.push(max_step);
            }
        }
        let points: Vec<Option<usize>> = beam_on
            .iter()
            .map(|&on| {
                if on {
                    o.points_per_chunk_on
                } else {
                    o.points_per_chunk_off
                }
            })
            .collect();
        Ok(make_chunks(
            &boundaries,
            &beam_on,
            &max_steps,
            beam,
            &points,
        ))
    }

    /// Run chunks one after the other, each from the previous final state.
    pub fn run_chunks(&mut self, chunks: &[SimChunk], o: &ChunkOptions) -> Result<OdeResult> {
        let mut res = self.new_result()?;
        for (i, chunk) in chunks.iter().enumerate() {
            let options = RunOptions {
                max_step: chunk.max_step_size,
                t_eval: None,
                ..o.run.clone()
            };
            let y0 = if i > 0 { res.final_cc() } else { o.y0.clone() };
            let sol = self.raw_run(
                (chunk.start, chunk.stop),
                &options,
                y0.as_deref(),
                chunk.beam_on,
                chunk.dose_rate,
                chunk.n_points,
            )?;
            let short = self.last_short_on_span;
            let mut dose_rate = if o.run.use_beam {
                chunk.dose_rate
            } else {
                None
            };
            if short.is_some() {
                dose_rate = Some(0.0);
            }
            let thinned = chunk.n_points.is_some() && short.is_none();
            res.append(
                sol,
                self.last_exec_time,
                &self.last_bursts(),
                dose_rate,
                short,
                thinned,
            );
        }
        Ok(res)
    }

    /// `prepare_chunked_run` then `run_chunks`.
    pub fn chunked_run(&mut self, t_span: (f64, f64), o: &ChunkOptions) -> Result<OdeResult> {
        let chunks = self.prepare_chunked_run(t_span, o)?;
        self.run_chunks(&chunks, o)
    }

    /// `prepare_single_pulse` then `run_chunks`.
    pub fn single_pulse_run(&mut self, t_span: (f64, f64), o: &ChunkOptions) -> Result<OdeResult> {
        let chunks = self.prepare_single_pulse(t_span, o)?;
        self.run_chunks(&chunks, o)
    }
}
