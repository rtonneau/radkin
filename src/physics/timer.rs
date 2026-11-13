//! Beam timers (port of radiopyo `physics/timer.py`).

use crate::consts::COUNT_RTOL;
use crate::error::{RadkinError, Result};
use crate::pyfloat::{floordiv, pymod};

/// State of a timer at a given time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimerState {
    Off,
    On,
}

/// `(n, remainder)` with `duration = n * period + remainder`, tolerant to
/// rounding at period boundaries.
pub fn split_periods(duration: f64, period: f64) -> (f64, f64) {
    let n = (duration / period + COUNT_RTOL).floor();
    let mut remainder = duration - n * period;
    if remainder <= COUNT_RTOL * period {
        remainder = 0.0;
    }
    (n, remainder)
}

/// Periodic ON/OFF timer; `period = on_time = inf` is a continuous beam.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Timer {
    pub period: f64,
    pub on_time: f64,
    pub start_time: f64,
}

impl Timer {
    /// Continuous timer.
    pub fn constant(start_time: Option<f64>) -> Self {
        Self {
            period: f64::INFINITY,
            on_time: f64::INFINITY,
            start_time: start_time.unwrap_or(0.0),
        }
    }

    /// Pulsed timer.
    pub fn pulsed(period: f64, on_time: f64, start_time: Option<f64>) -> Self {
        Self {
            period,
            on_time,
            start_time: start_time.unwrap_or(0.0),
        }
    }

    pub fn period_index(&self, time: f64) -> f64 {
        floordiv(time - self.start_time, self.period).trunc()
    }

    fn phase(&self, time: f64) -> f64 {
        time - self.start_time - self.period * self.period_index(time)
    }

    pub fn time_till_next_period(&self, time: f64) -> f64 {
        self.period - self.phase(time)
    }

    pub fn state_at(&self, time: f64) -> TimerState {
        if time >= self.start_time && pymod(time - self.start_time, self.period) <= self.on_time {
            TimerState::On
        } else {
            TimerState::Off
        }
    }

    pub fn time_till_next_on(&self, time: f64) -> f64 {
        if self.state_at(time) == TimerState::On {
            return 0.0;
        }
        self.time_till_next_period(time)
    }

    pub fn time_till_next_off(&self, time: f64) -> f64 {
        if self.state_at(time) == TimerState::Off {
            return 0.0;
        }
        self.on_time - self.phase(time)
    }

    pub fn is_continuous(&self) -> bool {
        (self.period == f64::INFINITY && self.on_time == f64::INFINITY)
            || self.on_time >= self.period
    }

    pub fn duty_cycle(&self) -> f64 {
        if self.is_continuous() {
            1.0
        } else {
            self.on_time / self.period
        }
    }
}

/// Trains of micro pulses: a macro timer whose ON phases run a sub-timer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MacroTimer {
    pub period: f64,
    pub on_time: f64,
    pub start_time: f64,
    pub sub: Timer,
}

impl MacroTimer {
    /// radiopyo `MacroTimer.new_macropulse`.
    pub fn new(
        macro_period: f64,
        macro_on_time: f64,
        micro_period: f64,
        micro_on_time: f64,
        macro_start_time: Option<f64>,
        micro_start_time: Option<f64>,
    ) -> Result<Self> {
        if macro_period <= micro_period {
            return Err(RadkinError::Beam(
                "Macro Period is <= than micro period. Do you really need a macro timer?".into(),
            ));
        }
        if let Some(ms) = micro_start_time
            && ms >= macro_on_time
        {
            return Err(RadkinError::Beam(format!(
                "micro_start_time ({ms}) must be < the train on_time ({macro_on_time}): no micro pulse would fit in a train."
            )));
        }
        Ok(Self {
            period: macro_period,
            on_time: macro_on_time,
            start_time: macro_start_time.unwrap_or(0.0),
            sub: Timer::pulsed(micro_period, micro_on_time, micro_start_time),
        })
    }

    fn period_index(&self, time: f64) -> f64 {
        floordiv(time - self.start_time, self.period).trunc()
    }

    fn phase(&self, time: f64) -> f64 {
        time - self.start_time - self.period * self.period_index(time)
    }

    /// Time relative to the current train start (`time >= start_time`).
    fn relative(&self, time: f64) -> f64 {
        self.phase(time)
    }

    pub fn state_at(&self, time: f64) -> TimerState {
        if time < self.start_time {
            return TimerState::Off;
        }
        let rel = self.relative(time);
        if rel <= self.on_time {
            self.sub.state_at(rel)
        } else {
            TimerState::Off
        }
    }

    pub fn micro_pulses_per_train(&self) -> f64 {
        let available = self.on_time - self.sub.start_time;
        if available <= 0.0 {
            return 0.0;
        }
        let (n, remainder) = split_periods(available, self.sub.period);
        let mut out = n;
        if remainder > 0.0 {
            out += (remainder / self.sub.on_time).min(1.0);
        }
        out
    }

    pub fn duty_cycle(&self) -> f64 {
        self.micro_pulses_per_train() * self.sub.on_time / self.period
    }

    pub fn time_till_next_on(&self, time: f64) -> f64 {
        if time < self.start_time {
            return self.start_time - time + self.sub.start_time;
        }
        if self.state_at(time) == TimerState::On {
            return 0.0;
        }
        let rel = self.relative(time);
        if rel <= self.on_time {
            let wait = self.sub.time_till_next_on(rel);
            if rel + wait < self.on_time {
                return wait;
            }
        }
        self.period - rel + self.sub.start_time
    }

    pub fn time_till_next_off(&self, time: f64) -> f64 {
        if self.state_at(time) == TimerState::Off {
            return 0.0;
        }
        let rel = self.relative(time);
        self.sub.time_till_next_off(rel).min(self.on_time - rel)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pulsed_timer_states() {
        let t = Timer::pulsed(1.0, 0.25, Some(0.5));
        assert_eq!(t.state_at(0.4), TimerState::Off);
        assert_eq!(t.state_at(0.6), TimerState::On);
        assert_eq!(t.state_at(0.8), TimerState::Off);
        assert_eq!(t.time_till_next_on(0.8), 0.7);
        assert!((t.time_till_next_off(0.6) - 0.15).abs() < 1e-15);
        assert_eq!(Timer::constant(None).state_at(1e9), TimerState::On);
    }

    #[test]
    fn macro_timer_counts_pulses() {
        let m = MacroTimer::new(0.1, 0.001, 1e-4, 2e-6, Some(0.01), Some(0.0)).unwrap();
        assert_eq!(m.micro_pulses_per_train(), 10.0);
        assert_eq!(m.state_at(0.0100001), TimerState::On);
        assert_eq!(m.state_at(0.0100500), TimerState::Off);
        assert!(MacroTimer::new(0.1, 0.001, 0.2, 2e-6, None, None).is_err());
    }
}
