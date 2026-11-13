use radkin::RadkinError;
use radkin::physics::Beam;
use radkin::unit_cell::{ChunkOptions, RunOptions, UnitCell, thinned_times};
use serde_json::Value;

fn runs() -> Value {
    serde_json::from_str(&std::fs::read_to_string("tests/golden/runs.json").unwrap()).unwrap()
}

fn pure_water() -> UnitCell {
    UnitCell::from_toml("data/advanced_water_radiolysis.toml").unwrap()
}

#[test]
fn thinned_times_match_radiopyo() {
    for case in runs()["thinned_times"].as_array().unwrap() {
        let a: Vec<f64> = case["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        let got = thinned_times(a[0], a[1], a[2] as usize);
        let exp: Vec<f64> = case["times"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        assert_eq!(got.len(), exp.len(), "{a:?}");
        for (g, e) in got.iter().zip(&exp) {
            assert!(
                (g - e).abs() <= 1e-15 * e.abs(),
                "{a:?}: {got:?} vs {exp:?}"
            );
        }
    }
}

#[test]
fn chunks_follow_the_golden_chunked_run() {
    let mut cell = pure_water();
    let beam = Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None).unwrap();
    cell.set_beam(beam, "default", Some(2.5)).unwrap();
    let opts = ChunkOptions {
        max_step_size_on: Some(1e-7),
        max_step_size_off: Some(1e-3),
        ..Default::default()
    };
    let chunks = cell.prepare_chunked_run((0.0, 0.05), &opts).unwrap();
    let g = &runs()["chunked_run"];
    assert_eq!(chunks.len(), g["n_points"].as_array().unwrap().len());
    let times: Vec<f64> = g["time"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    for c in &chunks {
        assert!(
            times.iter().any(|t| (t - c.stop).abs() <= 1e-12 * c.stop),
            "chunk stop {} not in golden",
            c.stop
        );
        assert_eq!(c.max_step_size, Some(if c.beam_on { 1e-7 } else { 1e-3 }));
    }
}

#[test]
fn span_shorter_than_min_on_time_is_added_analytically() {
    let mut cell = pure_water();
    cell.set_beam(
        Beam::constant(1.0, None, None).unwrap(),
        "default",
        Some(0.201),
    )
    .unwrap();
    let res = cell
        .run((1e-7, 1.5e-7), &RunOptions::default(), None)
        .unwrap();
    assert_eq!(res.t(), &[1.5e-7]);
    assert_eq!(res.short_on_spans.len(), 1);
    let ie = cell.env.species.index_of("e_aq").unwrap();
    let kr = cell
        .env
        .reactions
        .radiolytic
        .iter()
        .find(|r| r.species == "e_aq")
        .unwrap()
        .kr;
    let expected = kr * 1.0 * 5e-8 * 1e6;
    assert!((res.y()[0][ie] - expected).abs() <= 1e-12 * expected);
}

#[test]
fn t_eval_before_the_burst_is_dropped() {
    let mut cell = pure_water();
    cell.set_beam(
        Beam::constant(1.0, None, None).unwrap(),
        "default",
        Some(0.201),
    )
    .unwrap();
    let opts = RunOptions {
        t_eval: Some(vec![1e-7, 1.5e-7, 2e-4, 1e-3]),
        ..Default::default()
    };
    let res = cell.run((1e-7, 1e-3), &opts, None).unwrap();
    assert_eq!(res.t(), &[2e-4, 1e-3]);
}

#[test]
fn last_pulse_cut_by_max_dose_is_a_short_on_span() {
    let mut cell = pure_water();
    let beam = Beam::pulsed(1e6, 0.01, 2e-6, Some(4.3), None).unwrap();
    cell.set_beam(beam, "default", Some(0.201)).unwrap();
    let res = cell
        .chunked_run((0.0, 0.03), &ChunkOptions::default())
        .unwrap();
    assert_eq!(res.short_on_spans.len(), 1);
    let s = &res.short_on_spans[0];
    assert!(
        (s.start - 0.02).abs() < 1e-12 && (s.stop - s.start - 3e-7).abs() < 1e-12,
        "{s:?}"
    );
    assert_eq!(res.burst_dose_rates.as_deref(), Some(&[1e6, 1e6][..]));
}

#[test]
fn no_beam_is_an_error() {
    let cell = pure_water();
    let err = cell
        .clone()
        .run((0.0, 1.0), &RunOptions::default(), None)
        .unwrap_err();
    assert!(matches!(err, RadkinError::NoBeamDefined), "{err}");
}

#[test]
fn from_config_with_gvalues_matches_from_config() {
    let cfg = radkin::config::Config::load("data/config_Labarbe.toml").unwrap();
    let gv = radkin::gvalues::GValuesBoscolo::load(&radkin::consts::data_dir()).unwrap();
    let a = UnitCell::from_config(&cfg).unwrap();
    let b = UnitCell::from_config_with_gvalues(&cfg, gv).unwrap();
    assert_eq!(a.env, b.env);
    assert_eq!(a.mat, b.mat);
}
