//! The shortest beam ON time is a process-wide setting: these checks change it, so they
//! live in their own test binary and run in one test, in order.

use radkin::RadkinError;
use radkin::consts::{MIN_BEAM_ON_TIME, min_beam_on_time, set_min_beam_on_time};
use radkin::physics::Beam;
use radkin::unit_cell::{ChunkOptions, UnitCell};

fn last_pulse_run(max_dose: f64) -> usize {
    let mut cell = UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap();
    let beam = Beam::pulsed(1e6, 0.01, 2e-6, Some(max_dose), None).unwrap();
    cell.set_beam(beam, "default", Some(0.201)).unwrap();
    let res = cell
        .chunked_run((0.0, 0.03), &ChunkOptions::default())
        .unwrap();
    res.short_on_spans.len()
}

#[test]
fn min_beam_on_time_setting() {
    assert_eq!(min_beam_on_time(), MIN_BEAM_ON_TIME);
    for bad in [0.0, -1e-6, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            set_min_beam_on_time(Some(bad)),
            Err(RadkinError::Value(_))
        ));
    }
    assert_eq!(min_beam_on_time(), MIN_BEAM_ON_TIME);

    // default: a 0.3 us pulse is rejected, a 0.3 us last pulse is analytic
    assert!(matches!(
        Beam::pulsed(1e6, 0.01, 3e-7, None, None),
        Err(RadkinError::ShortBeamOn(_))
    ));
    assert_eq!(last_pulse_run(4.3), 1);

    // lower threshold: both are now accepted / integrated
    set_min_beam_on_time(Some(2e-7)).unwrap();
    assert_eq!(min_beam_on_time(), 2e-7);
    Beam::pulsed(1e6, 0.01, 3e-7, None, None).unwrap();
    assert_eq!(last_pulse_run(4.3), 0);

    // below T_ZERO: a span of T_ZERO (0.1 us last pulse) is still analytic
    set_min_beam_on_time(Some(1e-8)).unwrap();
    assert_eq!(last_pulse_run(4.1), 1);

    // higher threshold: rejects a 1 us pulse, names the setting
    set_min_beam_on_time(Some(1e-5)).unwrap();
    let err = Beam::pulsed(1e6, 0.01, 2e-6, None, None).unwrap_err();
    assert!(
        err.to_string().contains("MIN_BEAM_ON_TIME=1.000e-5"),
        "{err}"
    );

    set_min_beam_on_time(None).unwrap();
    assert_eq!(min_beam_on_time(), MIN_BEAM_ON_TIME);
}
