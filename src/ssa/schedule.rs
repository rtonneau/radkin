//! Dose-rate schedule of a run: the beam as constant-dose-rate segments.

use crate::error::{RadkinError, Result};
use crate::physics::Beam;
use crate::ssa::exact::Segment;

/// Segments of `t_span` at a constant dose rate, cut at the beam
/// transitions; one segment without dose when `beam` is `None`.
///
/// Each G-value then acts as a zero-order source reaction of propensity
/// `yield * dose_rate * N_A V`, constant over a segment.
pub fn beam_segments(beam: Option<&Beam>, t_span: (f64, f64)) -> Result<Vec<Segment>> {
    let (start, stop) = t_span;
    if !(start.is_finite() && stop.is_finite() && stop > start) {
        return Err(RadkinError::Value(format!(
            "ssa: t_span must be finite with stop > start, got ({start}, {stop})"
        )));
    }
    let Some(beam) = beam else {
        return Ok(vec![Segment {
            start,
            stop,
            dose_rate: 0.0,
        }]);
    };
    // Rows `[time, dose rate from that time on]`; the last row is `stop`.
    let rows = beam.make_matrix(stop, start);
    let mut segments: Vec<Segment> = Vec::with_capacity(rows.len());
    for w in rows.windows(2) {
        let (a, b) = (w[0][0], w[1][0]);
        if b <= a {
            continue;
        }
        match segments.last_mut() {
            // Merge neighbours at the same dose rate.
            Some(last) if last.dose_rate == w[0][1] => last.stop = b,
            _ => segments.push(Segment {
                start: a,
                stop: b,
                dose_rate: w[0][1],
            }),
        }
    }
    match (segments.first(), segments.last()) {
        (Some(f), Some(l)) if f.start == start && l.stop == stop => Ok(segments),
        _ => Err(RadkinError::Value(format!(
            "ssa: the beam schedule does not cover ({start}, {stop})"
        ))),
    }
}

/// Dose [Gy] delivered over `segments`.
pub fn delivered_dose(segments: &[Segment]) -> f64 {
    segments
        .iter()
        .map(|s| s.dose_rate * (s.stop - s.start))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_beam_is_one_segment_without_dose() {
        let s = beam_segments(None, (0.0, 2.0)).unwrap();
        assert_eq!(s.len(), 1);
        assert_eq!((s[0].start, s[0].stop, s[0].dose_rate), (0.0, 2.0, 0.0));
        assert!(beam_segments(None, (1.0, 1.0)).is_err());
        assert!(beam_segments(None, (0.0, f64::INFINITY)).is_err());
    }

    #[test]
    fn constant_beam_with_a_late_start() {
        let beam = Beam::constant(2.0, None, Some(1.0)).unwrap();
        let s = beam_segments(Some(&beam), (0.0, 3.0)).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!((s[0].stop, s[0].dose_rate), (1.0, 0.0));
        assert_eq!((s[1].stop, s[1].dose_rate), (3.0, 2.0));
        assert_eq!(delivered_dose(&s), 4.0);
    }

    #[test]
    fn pulsed_beam_alternates_and_delivers_its_dose() {
        // 3 pulses of 2 Gy: 1e6 Gy/s for 2 µs every 10 ms.
        let beam =
            Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
        let s = beam_segments(Some(&beam), (0.0, 0.05)).unwrap();
        let on: Vec<&Segment> = s.iter().filter(|x| x.dose_rate > 0.0).collect();
        assert_eq!(on.len(), 3);
        for (i, x) in on.iter().enumerate() {
            assert!((x.start - 0.01 * i as f64).abs() < 1e-12);
            assert!((x.stop - x.start - 2e-6).abs() < 1e-12);
            assert_eq!(x.dose_rate, 1e6);
        }
        assert!((delivered_dose(&s) - 6.0).abs() < 1e-6);
        assert!(s.windows(2).all(|w| w[0].stop == w[1].start));
        assert_eq!(s.last().unwrap().stop, 0.05);
    }

    #[test]
    fn span_starting_inside_a_pulse() {
        let beam =
            Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
        let s = beam_segments(Some(&beam), (1e-6, 5e-3)).unwrap();
        assert_eq!(s[0].dose_rate, 1e6);
        assert!((s[0].stop - 2e-6).abs() < 1e-15);
        assert!((delivered_dose(&s) - 1.0).abs() < 1e-9);
    }
}
