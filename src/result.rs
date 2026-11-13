//! Results of a run (port of radiopyo `simulation/ode_result.py`).

use crate::chemistry::env::SimEnv;
use crate::matrices::SimMatrices;
use crate::physics::Beam;
use crate::solver::Solution;

/// A beam ON span too short to integrate: its yield was added analytically.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct ShortOnSpan {
    pub start: f64,
    pub stop: f64,
    pub dose_rate: f64,
}

/// Concentrations over time plus what is needed to rebuild reaction rates.
#[derive(Debug, Clone)]
pub struct OdeResult {
    pub env: SimEnv,
    pub mat: SimMatrices,
    pub beam: Beam,
    pub let_: Option<f64>,
    pub t_zero_pulse: f64,
    t: Vec<f64>,
    y: Vec<Vec<f64>>,
    /// Dose rates of the analytic initial bursts (`None`: never set).
    pub burst_dose_rates: Option<Vec<f64>>,
    /// Constant dose rate of each appended block (`None`: the beam).
    pub step_dose_rates: Vec<Option<f64>>,
    pub short_on_spans: Vec<ShortOnSpan>,
    pub thinned: bool,
    pub n_points: Vec<usize>,
    pub nfev: Vec<usize>,
    pub njev: Vec<usize>,
    pub nst: Vec<usize>,
    pub exec_time: Vec<f64>,
    #[allow(dead_code)]
    pub(crate) thinned_warned: std::cell::Cell<bool>,
}

impl OdeResult {
    pub fn new(
        env: SimEnv,
        mat: SimMatrices,
        beam: Beam,
        let_: Option<f64>,
        t_zero_pulse: f64,
    ) -> Self {
        Self {
            env,
            mat,
            beam,
            let_,
            t_zero_pulse,
            t: Vec::new(),
            y: Vec::new(),
            burst_dose_rates: None,
            step_dose_rates: Vec::new(),
            short_on_spans: Vec::new(),
            thinned: false,
            n_points: Vec::new(),
            nfev: Vec::new(),
            njev: Vec::new(),
            nst: Vec::new(),
            exec_time: Vec::new(),
            thinned_warned: std::cell::Cell::new(false),
        }
    }

    /// Append one solved span (radiopyo `append_sol`).
    #[allow(clippy::too_many_arguments)]
    pub fn append(
        &mut self,
        sol: Solution,
        exec_time: f64,
        bursts: &[f64],
        dose_rate: Option<f64>,
        short: Option<ShortOnSpan>,
        thinned: bool,
    ) {
        self.thinned |= thinned;
        if let Some(s) = short {
            self.short_on_spans.push(s);
        }
        self.burst_dose_rates
            .get_or_insert_with(Vec::new)
            .extend_from_slice(bursts);
        self.n_points.push(sol.t.len());
        self.nfev.push(sol.nfev);
        self.njev.push(sol.njev);
        self.nst.push(sol.nst);
        self.t.extend(sol.t);
        self.y.extend(sol.y);
        self.step_dose_rates.push(dose_rate);
        self.exec_time.push(exec_time);
    }

    /// Stored times [s].
    pub fn t(&self) -> &[f64] {
        &self.t
    }

    /// Dynamic state [µmol/L] at each stored time.
    pub fn y(&self) -> &[Vec<f64>] {
        &self.y
    }

    /// Last stored state.
    pub fn final_cc(&self) -> Option<Vec<f64>> {
        self.y.last().cloned()
    }
}
