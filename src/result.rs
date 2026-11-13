//! Results of a run (port of radiopyo `simulation/ode_result.py`).

use crate::chemistry::env::SimEnv;
use crate::error::{RadkinError, Result};
use crate::kernel::{ReactionRates, Workspace, reaction_rates, to_full_molar};
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

/// numpy `trapezoid(y, x)`.
pub fn trapezoid(y: &[f64], x: &[f64]) -> f64 {
    (1..x.len())
        .map(|i| (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2.0)
        .sum()
}

/// Descending stable sort by value (radiopyo `sort_values(ascending=False)`).
fn sorted_desc(mut v: Vec<(String, f64)>) -> Vec<(String, f64)> {
    v.sort_by(|a, b| b.1.total_cmp(&a.1));
    v
}

impl OdeResult {
    fn warn_if_thinned(&self) {
        if self.thinned && !self.thinned_warned.get() {
            self.thinned_warned.set(true);
            log::warn!(
                "Integral over the stored points of a thinned result (points_per_chunk_on / points_per_chunk_off): coarser than over the solver steps; use more points per chunk, or no thinning, for accurate integrals"
            );
        }
    }

    /// Species labels (dynamic sorted, then fixed sorted) and rows of
    /// concentrations [µmol/L] (radiopyo `to_pandas`).
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
            .y
            .iter()
            .map(|y| {
                let mut r = y.clone();
                r.extend(&cst);
                r
            })
            .collect();
        (labels, rows)
    }

    /// Dose rate used for each stored point by the rate reconstruction.
    fn point_dose_rates(&self) -> Vec<f64> {
        let consistent = self.step_dose_rates.len() == self.n_points.len()
            && self.n_points.iter().sum::<usize>() == self.t.len();
        if !consistent || self.step_dose_rates.iter().all(Option::is_none) {
            return self.t.iter().map(|&t| self.beam.at(t)).collect();
        }
        let mut out = Vec::with_capacity(self.t.len());
        let mut first = 0;
        for (n, dr) in self.n_points.iter().zip(&self.step_dose_rates) {
            for &t in &self.t[first..first + n] {
                out.push(dr.unwrap_or_else(|| self.beam.at(t)));
            }
            first += n;
        }
        out
    }

    /// Reaction labels (k with pairs first, radiolytic, enzymatic) and rows
    /// of rates [mol/L/s] (radiopyo `reconstruct_reactions`).
    pub fn reaction_rates(&self) -> (Vec<String>, Vec<Vec<f64>>) {
        let mut ws = Workspace::new(&self.mat);
        let mut rates = ReactionRates::new(&self.mat);
        let mut y_full = vec![0.0; self.mat.n_species];
        let rows = self
            .y
            .iter()
            .zip(self.point_dose_rates())
            .map(|(y, dr)| {
                to_full_molar(y, &self.mat, &mut y_full);
                reaction_rates(&y_full, dr, &self.mat, &mut ws, &mut rates);
                rates.stacked()
            })
            .collect();
        (self.env.reaction_labels(), rows)
    }

    /// Indices of the stored points with `start <= t <= stop`.
    fn window(&self, start: Option<f64>, stop: Option<f64>) -> Result<Vec<usize>> {
        let (a, b) = (
            start.unwrap_or(f64::NEG_INFINITY),
            stop.unwrap_or(f64::INFINITY),
        );
        if a >= b {
            return Err(RadkinError::Value(format!(
                "Reversed boundaries. Trying to slice from {a} s to {b} s"
            )));
        }
        Ok((0..self.t.len())
            .filter(|&i| self.t[i] >= a && self.t[i] <= b)
            .collect())
    }

    /// Time integral of every species over `start <= t < stop` [µmol/L·s].
    pub fn integrate_species(&self, start: f64, stop: f64) -> Vec<(String, f64)> {
        self.warn_if_thinned();
        let (labels, rows) = self.concentrations();
        let idx: Vec<usize> = (0..self.t.len())
            .filter(|&i| self.t[i] >= start && self.t[i] < stop)
            .collect();
        let x: Vec<f64> = idx.iter().map(|&i| self.t[i]).collect();
        labels
            .into_iter()
            .enumerate()
            .map(|(j, l)| {
                let y: Vec<f64> = idx.iter().map(|&i| rows[i][j]).collect();
                (l, trapezoid(&y, &x))
            })
            .collect()
    }

    fn burst_list(&self) -> Vec<f64> {
        match &self.burst_dose_rates {
            Some(b) => b.clone(),
            None => vec![self.beam.tot_pulses() * self.beam.dose_rate_peak],
        }
    }

    /// Net yield series of the pairs involving `species` (radiopyo
    /// `compute_yields_of_eq_reactions_by_species`).
    fn eq_yields(
        &self,
        species: &str,
        rates: &[Vec<f64>],
        idx: &[usize],
    ) -> Vec<(String, Vec<f64>)> {
        let mut out = Vec::new();
        for bi in self.env.reactions.iter_pairs() {
            if !(bi.forward.contains(species) || bi.backward.contains(species)) {
                continue;
            }
            let (fwd, bwd) = if bi.backward.has_product(species) {
                (&bi.backward, &bi.forward)
            } else {
                (&bi.forward, &bi.backward)
            };
            let sf = fwd.stoi_of(species).unwrap_or(0) as f64;
            let sb = bwd.stoi_of(species).unwrap_or(0) as f64;
            let series = idx
                .iter()
                .map(|&i| sf * rates[i][fwd.index] - sb * rates[i][bwd.index])
                .collect();
            out.push((fwd.as_label(), series));
        }
        out
    }

    fn column(labels: &[String], label: &str) -> usize {
        labels
            .iter()
            .position(|l| l == label)
            .expect("reaction label")
    }

    fn integrate_columns(cols: Vec<(String, Vec<f64>)>, x: &[f64]) -> Vec<(String, f64)> {
        cols.into_iter()
            .map(|(l, s)| {
                let y: Vec<f64> = s.iter().map(|v| v * 1e6).collect();
                (l, trapezoid(&y, x))
            })
            .collect()
    }

    /// Gross production of `species` per reaction [µmol/L], sorted
    /// descending (radiopyo `integrate_production_of`).
    pub fn integrate_production_of(
        &self,
        species: &str,
        start: Option<f64>,
        stop: Option<f64>,
        include_equilibrium: bool,
        include_radiolytic: bool,
    ) -> Result<Vec<(String, f64)>> {
        self.warn_if_thinned();
        let idx = self.window(start, stop)?;
        let (labels, rates) = self.reaction_rates();
        let x: Vec<f64> = idx.iter().map(|&i| self.t[i]).collect();
        let r = &self.env.reactions;
        let scaled = |label: String, stoi: u32| {
            let c = Self::column(&labels, &label);
            let s: Vec<f64> = idx.iter().map(|&i| rates[i][c] * stoi as f64).collect();
            (label, s)
        };
        let mut cols: Vec<(String, Vec<f64>)> = Vec::new();
        for k in r.k.iter().filter(|k| k.has_product(species)) {
            cols.push(scaled(k.as_label(), k.stoi_of(species).unwrap_or(0)));
        }
        for e in r.enzymatic.iter().filter(|e| e.has_product(species)) {
            cols.push(scaled(e.as_label(), e.stoi_of(species).unwrap_or(0)));
        }
        if include_equilibrium {
            for (l, s) in self.eq_yields(species, &rates, &idx) {
                if trapezoid(&s, &x) >= 0.0 {
                    cols.push((l, s));
                }
            }
        }
        let radiolytic: Vec<_> = r
            .radiolytic
            .iter()
            .filter(|g| g.species == species)
            .collect();
        if include_radiolytic {
            for g in &radiolytic {
                cols.push(scaled(g.label(), 1));
            }
        }
        let mut out = Self::integrate_columns(cols, &x);
        if include_radiolytic {
            for g in radiolytic {
                let bursts: f64 = self
                    .burst_list()
                    .iter()
                    .map(|dr| g.initial_cc_before(self.t_zero_pulse, *dr) * 1e6)
                    .sum();
                let spans: f64 = self
                    .short_on_spans
                    .iter()
                    .map(|s| g.initial_cc_before(s.stop - s.start, s.dose_rate) * 1e6)
                    .sum();
                if let Some(e) = out.iter_mut().find(|(l, _)| *l == g.label()) {
                    e.1 = e.1 + bursts + spans;
                }
            }
        }
        Ok(sorted_desc(out))
    }

    /// Gross consumption of `species` per reaction [µmol/L] (negative
    /// values), sorted descending (radiopyo `integrate_consumption_of`).
    pub fn integrate_consumption_of(
        &self,
        species: &str,
        start: Option<f64>,
        stop: Option<f64>,
        include_equilibrium: bool,
    ) -> Result<Vec<(String, f64)>> {
        self.warn_if_thinned();
        let idx = self.window(start, stop)?;
        let (labels, rates) = self.reaction_rates();
        let x: Vec<f64> = idx.iter().map(|&i| self.t[i]).collect();
        let r = &self.env.reactions;
        let scaled = |label: String, stoi: u32| {
            let c = Self::column(&labels, &label);
            let s: Vec<f64> = idx.iter().map(|&i| -rates[i][c] * stoi as f64).collect();
            (label, s)
        };
        let mut cols: Vec<(String, Vec<f64>)> = Vec::new();
        for k in r.k.iter().filter(|k| k.has_reactant(species)) {
            cols.push(scaled(k.as_label(), k.stoi_of(species).unwrap_or(0)));
        }
        for e in r.enzymatic.iter().filter(|e| e.has_reactant(species)) {
            cols.push(scaled(e.as_label(), e.stoi_of(species).unwrap_or(0)));
        }
        if include_equilibrium {
            for (l, s) in self.eq_yields(species, &rates, &idx) {
                if trapezoid(&s, &x) <= 0.0 {
                    cols.push((l, s));
                }
            }
        }
        Ok(sorted_desc(Self::integrate_columns(cols, &x)))
    }
}
