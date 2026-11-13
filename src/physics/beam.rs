//! Particle beams (port of radiopyo `physics/beam.py`).

use crate::consts::{COUNT_RTOL, MIN_BEAM_ON_TIME, ON_TIME_RTOL};
use crate::error::{RadkinError, Result};
use crate::physics::timer::{MacroTimer, Timer, TimerState, split_periods};
use crate::pyfloat::{floordiv, pymod, round6};

/// Beam type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeamKind {
    Constant,
    SinglePulse,
    Pulsed,
    MacroPulsed(MacroTimer),
}

/// A beam: peak dose rate, time structure and maximum dose.
#[derive(Debug, Clone, PartialEq)]
pub struct Beam {
    pub kind: BeamKind,
    /// Dose rate while ON [Gy/s].
    pub dose_rate_peak: f64,
    /// Timer (the micro-pulse timer of a macro-pulsed beam).
    pub timer: Timer,
    /// Maximum dose [Gy], `inf` if none.
    pub max_dose: f64,
    /// Time at which `max_dose` is reached [s].
    pub time_at_max_dose: f64,
    /// Average dose rate [Gy/s].
    pub dose_rate: f64,
}

fn isclose(a: f64, b: f64, rel: f64) -> bool {
    a == b || (a - b).abs() <= rel * a.abs().max(b.abs())
}

fn round2(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    format!("{x:.2}").parse().unwrap_or(x)
}

impl Beam {
    fn build(
        kind: BeamKind,
        dose_rate_peak: f64,
        timer: Timer,
        max_dose: Option<f64>,
    ) -> Result<Self> {
        let mut beam = Self {
            kind,
            dose_rate_peak,
            timer,
            max_dose: max_dose.unwrap_or(f64::INFINITY),
            time_at_max_dose: f64::INFINITY,
            dose_rate: dose_rate_peak,
        };
        beam.post_init()?;
        Ok(beam)
    }

    fn post_init(&mut self) -> Result<()> {
        if self.dose_rate_peak <= 0.0 {
            return Err(RadkinError::Beam(format!(
                "Null or Negative inst. dose rate for beam ({})",
                self.dose_rate_peak
            )));
        }
        self.time_at_max_dose = self.compute_time_at_max_dose(self.max_dose);
        self.dose_rate = self.compute_average_dose_rate();
        self.check_on_times()
    }

    /// Constant beam.
    pub fn constant(
        dose_rate: f64,
        max_dose: Option<f64>,
        start_time: Option<f64>,
    ) -> Result<Self> {
        Self::build(
            BeamKind::Constant,
            dose_rate,
            Timer::constant(start_time),
            max_dose,
        )
    }

    /// Single pulse delivering `max_dose` at `dose_rate`.
    pub fn single_pulse(dose_rate: f64, max_dose: f64, start_time: Option<f64>) -> Result<Self> {
        Self::build(
            BeamKind::SinglePulse,
            dose_rate,
            Timer::constant(start_time),
            Some(max_dose),
        )
    }

    /// Single pulse of duration `on_time`.
    pub fn single_pulse_from_on_time(
        dose_rate: f64,
        on_time: f64,
        start_time: Option<f64>,
    ) -> Result<Self> {
        Self::single_pulse(dose_rate, on_time * dose_rate, start_time)
    }

    /// Pulsed beam.
    pub fn pulsed(
        dose_rate_peak: f64,
        period: f64,
        on_time: f64,
        max_dose: Option<f64>,
        start_time: Option<f64>,
    ) -> Result<Self> {
        Self::build(
            BeamKind::Pulsed,
            dose_rate_peak,
            Timer::pulsed(period, on_time, start_time),
            max_dose,
        )
    }

    /// radiopyo `PulsedBeam.from_dose_rate_peak`.
    pub fn pulsed_from_dose_rate_peak(
        peak: f64,
        period: f64,
        on_time: f64,
        n_pulse: Option<f64>,
        max_dose: Option<f64>,
        start_time: Option<f64>,
    ) -> Result<Self> {
        let n_pulse = n_pulse.unwrap_or(f64::INFINITY);
        let max_dose_n_pulse = on_time * n_pulse * peak;
        let max_dose = max_dose.unwrap_or(f64::INFINITY);
        if max_dose.is_finite()
            && max_dose_n_pulse.is_finite()
            && !isclose(max_dose_n_pulse, max_dose, 1e-9)
        {
            return Err(RadkinError::Beam(format!(
                "Inconsistent max dose provided by both max dose and number of pulses (max dose: {max_dose} || DPP: {max_dose_n_pulse})"
            )));
        }
        Self::pulsed(
            peak,
            period,
            on_time,
            Some(max_dose.min(max_dose_n_pulse)),
            start_time,
        )
    }

    /// radiopyo `PulsedBeam.from_dose_per_pulse`.
    pub fn pulsed_from_dose_per_pulse(
        pulse_dose: f64,
        period: f64,
        on_time: f64,
        n_pulse: Option<f64>,
        max_dose: Option<f64>,
        start_time: Option<f64>,
    ) -> Result<Self> {
        Self::pulsed_from_dose_rate_peak(
            pulse_dose / on_time,
            period,
            on_time,
            n_pulse,
            max_dose,
            start_time,
        )
    }

    /// Macro-pulsed beam (trains of micro pulses).
    #[allow(clippy::too_many_arguments)]
    pub fn macro_pulsed(
        dose_rate_peak: f64,
        period: f64,
        on_time: f64,
        micro_period: f64,
        micro_on_time: f64,
        max_dose: Option<f64>,
        micro_start_time: Option<f64>,
        start_time: Option<f64>,
    ) -> Result<Self> {
        let mt = MacroTimer::new(
            period,
            on_time,
            micro_period,
            micro_on_time,
            start_time,
            micro_start_time,
        )?;
        Self::build(BeamKind::MacroPulsed(mt), dose_rate_peak, mt.sub, max_dose)
    }

    /// Lower-case type name, as in radiopyo `beam_type().lower()`.
    pub fn beam_type(&self) -> &'static str {
        match self.kind {
            BeamKind::Constant => "constant",
            BeamKind::SinglePulse => "singlepulse",
            BeamKind::Pulsed => "pulsed",
            BeamKind::MacroPulsed(_) => "macro_pulsed",
        }
    }

    pub fn dose_per_pulse(&self) -> f64 {
        self.dose_rate_peak * self.timer.on_time
    }

    fn base_time_at_max_dose(&self, max_dose: f64) -> f64 {
        let t = &self.timer;
        if max_dose >= f64::INFINITY || self.dose_rate_peak <= 0.0 {
            return f64::INFINITY;
        }
        if max_dose == 0.0 {
            return 0.0;
        }
        if t.period >= f64::INFINITY {
            return t.start_time + max_dose / self.dose_rate_peak;
        }
        let dpp = self.dose_per_pulse();
        let n = floordiv(max_dose, dpp);
        let last_frac = round6(pymod(max_dose, dpp) / dpp);
        if n <= 0.0 {
            return t.start_time + last_frac * t.on_time;
        }
        if last_frac <= 0.0 {
            return t.start_time + (n - 1.0) * t.period + t.on_time;
        }
        t.start_time + n * t.period + last_frac * t.on_time
    }

    fn compute_time_at_max_dose(&self, max_dose: f64) -> f64 {
        let BeamKind::MacroPulsed(m) = self.kind else {
            return self.base_time_at_max_dose(max_dose);
        };
        if max_dose == 0.0 {
            return 0.0;
        }
        if max_dose >= f64::INFINITY || self.dose_rate_peak <= 0.0 {
            return f64::INFINITY;
        }
        let dpt = self.dose_per_train();
        let n = (max_dose / dpt + COUNT_RTOL).floor();
        let remaining = max_dose - n * dpt;
        if remaining <= COUNT_RTOL * dpt {
            return m.start_time + (n - 1.0) * m.period + self.base_time_at_max_dose(dpt);
        }
        m.start_time + n * m.period + self.base_time_at_max_dose(remaining)
    }

    /// Dose of one train of a macro-pulsed beam.
    pub fn dose_per_train(&self) -> f64 {
        match self.kind {
            BeamKind::MacroPulsed(m) => self.dose_per_pulse() * m.micro_pulses_per_train(),
            _ => self.dose_per_pulse(),
        }
    }

    fn compute_average_dose_rate(&self) -> f64 {
        if let BeamKind::MacroPulsed(m) = self.kind {
            let duration = self.time_at_max_dose - m.start_time;
            if duration.is_infinite() || duration <= 0.0 {
                return self.dose_rate_peak * m.duty_cycle();
            }
            return self.dose_rate_peak * (self.tot_pulses() * self.timer.on_time / duration);
        }
        if self.timer.period.is_infinite() {
            return self.dose_rate_peak;
        }
        let duration = self.time_at_max_dose - self.timer.start_time;
        if duration.is_infinite() || duration <= 0.0 {
            return self.dose_rate_peak * self.timer.duty_cycle();
        }
        self.dose_rate_peak * (self.tot_pulses() * self.timer.on_time / duration)
    }

    /// Number of pulses delivered (fractional for a cut last pulse).
    pub fn tot_pulses(&self) -> f64 {
        if let BeamKind::MacroPulsed(_) = self.kind {
            return round2(self.max_dose / self.dose_per_pulse());
        }
        if self.timer.period.is_infinite() {
            return 1.0;
        }
        let duration = self.time_at_max_dose - self.timer.start_time;
        if duration.is_infinite() {
            return f64::INFINITY;
        }
        if duration <= 0.0 {
            return 0.0;
        }
        let (n, remainder) = split_periods(duration, self.timer.period);
        n + (remainder / self.timer.on_time).min(1.0)
    }

    fn check_on_times(&self) -> Result<()> {
        let threshold = MIN_BEAM_ON_TIME * (1.0 - ON_TIME_RTOL);
        let mut on_times: Vec<(&str, f64)> = Vec::new();
        match self.kind {
            BeamKind::MacroPulsed(m) => {
                on_times.push(("micro_on_time", self.timer.on_time));
                on_times.push(("on_time", m.on_time));
            }
            _ => {
                if self.timer.period < f64::INFINITY {
                    on_times.push(("on_time", self.timer.on_time));
                } else if self.max_dose > 0.0 && self.time_at_max_dose.is_finite() {
                    on_times.push(("on_time", self.time_at_max_dose - self.timer.start_time));
                }
            }
        }
        for (name, on_time) in on_times {
            if on_time < threshold {
                return Err(RadkinError::ShortBeamOn(format!(
                    "{} beam: beam ON time {name}={on_time:.3e} s is shorter than MIN_BEAM_ON_TIME={MIN_BEAM_ON_TIME:.3e} s",
                    self.beam_type()
                )));
            }
        }
        Ok(())
    }

    /// Dose rate at `time` [Gy/s].
    pub fn at(&self, time: f64) -> f64 {
        match self.kind {
            BeamKind::MacroPulsed(m) => {
                if time < self.time_at_max_dose && m.state_at(time) == TimerState::On {
                    self.dose_rate_peak
                } else {
                    0.0
                }
            }
            _ => {
                if self.timer.start_time <= time
                    && time < self.time_at_max_dose
                    && self.timer.state_at(time) == TimerState::On
                {
                    self.dose_rate_peak
                } else {
                    0.0
                }
            }
        }
    }

    fn first_on_start(&self) -> f64 {
        match self.kind {
            BeamKind::MacroPulsed(m) => m.start_time,
            _ => self.timer.start_time,
        }
    }

    /// ON/OFF transitions from `start` on (radiopyo `_transitions`).
    fn transitions(&self, start: f64) -> Box<dyn Iterator<Item = (f64, bool)>> {
        if let BeamKind::MacroPulsed(m) = self.kind {
            return Box::new(macro_transitions(m, start));
        }
        let timer = self.timer;
        if timer.period.is_infinite() {
            let mut v = Vec::new();
            if start < timer.start_time {
                v.push((start, false));
                v.push((timer.start_time, true));
            } else {
                v.push((start, true));
            }
            v.push((f64::INFINITY, false));
            return Box::new(v.into_iter());
        }
        let mut head = Vec::new();
        let mut t_on;
        if start < timer.start_time {
            head.push((start, false));
            t_on = timer.start_time;
        } else if timer.state_at(start) == TimerState::On {
            head.push((start, true));
            head.push((start + timer.time_till_next_off(start), false));
            t_on = start + timer.time_till_next_period(start);
        } else {
            head.push((start, false));
            t_on = start + timer.time_till_next_period(start);
        }
        let mut pending_off: Option<f64> = None;
        let tail = std::iter::from_fn(move || {
            if let Some(t) = pending_off.take() {
                return Some((t, false));
            }
            let t = t_on;
            pending_off = Some(t + timer.on_time);
            t_on = t + timer.period;
            Some((t, true))
        });
        Box::new(head.into_iter().chain(tail))
    }

    /// `[time, dose rate]` rows of the beam between `start` and `stop`.
    pub fn make_matrix(&self, stop: f64, start: f64) -> Vec<[f64; 2]> {
        if start > self.time_at_max_dose || stop < self.first_on_start() {
            return vec![[start, 0.0], [stop, 0.0]];
        }
        let end = stop.min(self.time_at_max_dose);
        let mut walk = self.transitions(start);
        let (start_time, start_on) = walk.next().expect("transitions start with `start`");
        let peak = self.dose_rate_peak;
        let mut rows = vec![[start_time, if start_on { peak } else { 0.0 }]];
        for (time, is_on) in walk {
            if time < end {
                rows.push([time, if is_on { peak } else { 0.0 }]);
            } else {
                if !is_on {
                    rows.push([end, 0.0]);
                }
                break;
            }
        }
        if rows.last().unwrap()[0] < stop {
            rows.push([stop, 0.0]);
        } else {
            rows.last_mut().unwrap()[1] = 0.0;
        }
        rows
    }

    /// As [`Beam::make_matrix`] with 1/0 for ON/OFF.
    pub fn make_time_intervals(&self, stop: f64, start: f64) -> Vec<[f64; 2]> {
        self.make_matrix(stop, start)
            .into_iter()
            .map(|[t, d]| [t, if d > 0.0 { 1.0 } else { 0.0 }])
            .collect()
    }

    /// Beam from a config table (radiopyo `from_dict` of each beam type).
    pub fn from_table(beam_type: &str, param: &toml::Table) -> Result<Self> {
        let get = |k: &str| -> Result<Option<f64>> {
            match param.get(k) {
                None => Ok(None),
                Some(toml::Value::Float(f)) => Ok(Some(*f)),
                Some(toml::Value::Integer(i)) => Ok(Some(*i as f64)),
                Some(v) => Err(RadkinError::Beam(format!(
                    "'{k}' must be a number, got {v}"
                ))),
            }
        };
        let has = |k: &str| param.contains_key(k);
        let need = |k: &str, msg: &str| -> Result<f64> {
            get(k)?.ok_or_else(|| RadkinError::Beam(msg.to_string()))
        };
        match beam_type.to_lowercase().as_str() {
            "constant" => {
                let dose_rate = if has("on_time") && has("max_dose") {
                    need("max_dose", "")? / need("on_time", "")?
                } else {
                    let mut dr = 0.0;
                    for key in ["dose_rate_peak", "dose_rate_avg", "dose_rate"] {
                        if let Some(v) = get(key)? {
                            dr = v;
                            break;
                        }
                    }
                    dr
                };
                Self::constant(dose_rate, get("max_dose")?, get("start_time")?)
            }
            "singlepulse" => {
                let dr = need(
                    "dose_rate",
                    "SinglePulseBeam must be defined with a 'dose_rate' parameter",
                )?;
                if let Some(md) = get("max_dose")? {
                    Self::single_pulse(dr, md, get("start_time")?)
                } else if let Some(on) = get("on_time")? {
                    Self::single_pulse_from_on_time(dr, on, get("start_time")?)
                } else {
                    Err(RadkinError::Beam(
                        "SinglePulseBeam must be defined with either one of 'on_time' or 'max_dose' parameter".into(),
                    ))
                }
            }
            "pulsed" => {
                let period = match (get("frequency")?, get("period")?) {
                    (Some(f), _) => 1.0 / f,
                    (None, Some(p)) => p,
                    _ => {
                        return Err(RadkinError::Beam(
                            "PulsedBeam must be defined with a 'period' or 'frequency' parameter"
                                .into(),
                        ));
                    }
                };
                let (pd, peak, on, dr) = (
                    get("pulse_dose")?,
                    get("dose_rate_peak")?,
                    get("on_time")?,
                    get("dose_rate")?,
                );
                let (peak, on, pulse_dose) = match (pd, peak, on, dr) {
                    (Some(pd), Some(peak), _, _) => (Some(peak), pd / peak, Some(pd)),
                    (Some(pd), None, Some(on), _) => (Some(pd / on), on, Some(pd)),
                    (None, Some(peak), Some(on), _) => (Some(peak), on, Some(peak * on)),
                    (None, None, Some(on), Some(dr)) => {
                        let timer = Timer::pulsed(period, on, None);
                        (Some(dr / timer.duty_cycle()), on, None)
                    }
                    _ => {
                        return Err(RadkinError::Beam(
                            "PulsedBeam must be defined with at least 2 parameters out of 'dose_rate_peak', 'on_time' or 'pulse_dose'".into(),
                        ));
                    }
                };
                let (n_pulse, max_dose, start) =
                    (get("n_pulse")?, get("max_dose")?, get("start_time")?);
                match (peak, pulse_dose) {
                    (Some(peak), _) => {
                        Self::pulsed_from_dose_rate_peak(peak, period, on, n_pulse, max_dose, start)
                    }
                    (None, Some(pd)) => {
                        Self::pulsed_from_dose_per_pulse(pd, period, on, n_pulse, max_dose, start)
                    }
                    _ => unreachable!(),
                }
            }
            "macro_pulsed" => {
                let period = match (get("frequency")?, get("period")?) {
                    (Some(f), _) => 1.0 / f,
                    (None, Some(p)) => p,
                    _ => {
                        return Err(RadkinError::Beam(
                            "MacroPulsedBeam must be defined with a 'period' or 'frequency' parameter".into(),
                        ));
                    }
                };
                let on = need(
                    "on_time",
                    "MacroPulsedBeam must be defined with an 'on_time' parameter",
                )?;
                let micro_period = match (get("micro_frequency")?, get("micro_period")?) {
                    (Some(f), _) => 1.0 / f,
                    (None, Some(p)) => p,
                    _ => {
                        return Err(RadkinError::Beam(
                            "MacroPulsedBeam must be defined with an 'micro_period' or 'micro_frequency' parameter".into(),
                        ));
                    }
                };
                let micro_on = need(
                    "micro_on_time",
                    "MacroPulsedBeam must be defined with an 'micro_on_time' parameter",
                )?;
                let peak = need(
                    "dose_rate_peak",
                    "MacroPulsedBeam must be defined with a 'dose_rate_peak' parameter",
                )?;
                Self::macro_pulsed(
                    peak,
                    period,
                    on,
                    micro_period,
                    micro_on,
                    get("max_dose")?,
                    get("micro_start_time")?,
                    get("start_time")?,
                )
            }
            other => Err(RadkinError::Beam(format!("Unknown beam type: {other}"))),
        }
    }

    /// Parameters of the beam (radiopyo `as_dict`).
    pub fn as_table(&self) -> toml::Table {
        let mut t = toml::Table::new();
        let mut put = |k: &str, v: f64| {
            t.insert(k.to_string(), toml::Value::Float(v));
        };
        match self.kind {
            BeamKind::Constant | BeamKind::SinglePulse => {
                put("dose_rate", self.dose_rate);
                put("max_dose", self.max_dose);
                put("start_time", self.timer.start_time);
            }
            BeamKind::Pulsed => {
                put("dose_rate_peak", self.dose_rate_peak);
                put("period", self.timer.period);
                put("on_time", self.timer.on_time);
                put("max_dose", self.max_dose);
                put("start_time", self.timer.start_time);
            }
            BeamKind::MacroPulsed(m) => {
                put("dose_rate_peak", self.dose_rate_peak);
                put("period", m.period);
                put("on_time", m.on_time);
                put("micro_period", self.timer.period);
                put("micro_on_time", self.timer.on_time);
                put("max_dose", self.max_dose);
                put("micro_start_time", self.timer.start_time);
                put("start_time", m.start_time);
            }
        }
        t
    }
}

/// radiopyo `MacroPulsedBeam._transitions`.
fn macro_transitions(m: MacroTimer, start: f64) -> impl Iterator<Item = (f64, bool)> {
    let micro = m.sub;
    let mut is_on = m.state_at(start) == TimerState::On;
    let first = (start, is_on);
    let mut train_idx = ((start - m.start_time) / m.period).floor().max(0.0);
    let tolerance = COUNT_RTOL * micro.on_time;
    let mut pulse_idx = 0.0;
    let mut queue: std::collections::VecDeque<(f64, bool)> = std::collections::VecDeque::new();
    std::iter::once(first).chain(std::iter::from_fn(move || {
        loop {
            if let Some(x) = queue.pop_front() {
                return Some(x);
            }
            let train_start = m.start_time + train_idx * m.period;
            let train_end = train_start + m.on_time;
            let t_on = train_start + micro.start_time + pulse_idx * micro.period;
            if t_on >= train_end - tolerance {
                train_idx += 1.0;
                pulse_idx = 0.0;
                continue;
            }
            let t_off = (t_on + micro.on_time).min(train_end);
            for (time, state) in [(t_on, true), (t_off, false)] {
                if time >= start && state != is_on {
                    is_on = state;
                    queue.push_back((time, state));
                }
            }
            pulse_idx += 1.0;
        }
    }))
}
