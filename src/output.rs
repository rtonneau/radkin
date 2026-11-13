//! CSV/JSON outputs mirroring radiopyo's DataFrames.
//!
//! | file | radiopyo |
//! |---|---|
//! | `concentrations.csv` | `res.to_pandas()` (µmol/L) |
//! | `reaction_rates.csv` | `res.reconstruct_reactions()` (mol/L/s) |
//! | `species_integrals.csv` | `res.integrate_species()` (µmol/L·s) |
//! | `production_<sp>.csv` | `res.integrate_production_of(sp)` (µmol/L) |
//! | `consumption_<sp>.csv` | `res.integrate_consumption_of(sp)` (µmol/L) |

use std::fmt::Write as _;
use std::path::Path;

use indexmap::IndexMap;
use serde::Serialize;

use crate::error::{RadkinError, Result};
use crate::result::{OdeResult, ShortOnSpan};

/// Revision of rsoda pinned in `Cargo.toml`.
pub const RSODA_REV: &str = "312b742";

/// Run metadata written to `run.json`.
#[derive(Debug, Clone, Serialize)]
pub struct RunMeta {
    pub config: String,
    pub mode: String,
    pub t_span: (f64, f64),
    pub beam_label: String,
    pub beam_type: String,
    pub beam_params: IndexMap<String, f64>,
    #[serde(rename = "LET")]
    pub let_: Option<f64>,
    pub temperature: f64,
    pub g_values: IndexMap<String, f64>,
    pub burst_dose_rates: Option<Vec<f64>>,
    pub short_on_spans: Vec<ShortOnSpan>,
    pub thinned: bool,
    pub n_points: Vec<usize>,
    pub nfev: Vec<usize>,
    pub njev: Vec<usize>,
    pub nst: Vec<usize>,
    pub exec_time: f64,
    pub radkin_version: String,
    pub rsoda_rev: String,
    pub species_labels: Vec<String>,
    /// `as_label() -> str(reaction)`.
    pub reaction_strings: IndexMap<String, String>,
}

impl RunMeta {
    /// Metadata of `res`.
    pub fn from_result(
        res: &OdeResult,
        config: &str,
        mode: &str,
        t_span: (f64, f64),
        beam_label: &str,
    ) -> Self {
        let beam_params = res
            .beam
            .as_table()
            .into_iter()
            .filter_map(|(k, v)| v.as_float().map(|f| (k, f)))
            .collect();
        let mut species_labels = res.env.species.sorted_dyn_labels();
        species_labels.extend(res.env.species.sorted_cst_labels());
        Self {
            config: config.to_string(),
            mode: mode.to_string(),
            t_span,
            beam_label: beam_label.to_string(),
            beam_type: res.beam.beam_type().to_string(),
            beam_params,
            let_: res.let_,
            temperature: res.env.temperature,
            g_values: res
                .env
                .reactions
                .radiolytic
                .iter()
                .map(|r| (r.label(), r.ge()))
                .collect(),
            burst_dose_rates: res.burst_dose_rates.clone(),
            short_on_spans: res.short_on_spans.clone(),
            thinned: res.thinned,
            n_points: res.n_points.clone(),
            nfev: res.nfev.clone(),
            njev: res.njev.clone(),
            nst: res.nst.clone(),
            exec_time: res.exec_time.iter().sum(),
            radkin_version: env!("CARGO_PKG_VERSION").to_string(),
            rsoda_rev: RSODA_REV.to_string(),
            species_labels,
            reaction_strings: res
                .env
                .reaction_labels()
                .into_iter()
                .zip(res.env.reaction_strings())
                .collect(),
        }
    }
}

fn write(path: &Path, text: &str) -> Result<()> {
    std::fs::write(path, text).map_err(|e| RadkinError::io(path, e))
}

/// Table with a `t` column then one column per label (full precision).
pub fn table_csv(labels: &[String], t: &[f64], rows: &[Vec<f64>]) -> String {
    let mut s = String::from("t");
    for l in labels {
        s.push(',');
        s.push_str(l);
    }
    s.push('\n');
    for (ti, row) in t.iter().zip(rows) {
        let _ = write!(s, "{ti:e}");
        for v in row {
            let _ = write!(s, ",{v:e}");
        }
        s.push('\n');
    }
    s
}

/// Two-column `key,value` series.
pub fn series_csv(key: &str, values: &[(String, f64)]) -> String {
    let mut s = format!("{key},value\n");
    for (k, v) in values {
        let _ = writeln!(s, "{k},{v:e}");
    }
    s
}

/// Write every output file of a run into `dir`.
pub fn write_run_dir(
    dir: &Path,
    res: &OdeResult,
    meta: &RunMeta,
    production: &[String],
    consumption: &[String],
) -> Result<()> {
    std::fs::create_dir_all(dir).map_err(|e| RadkinError::io(dir, e))?;
    let (labels, rows) = res.concentrations();
    write(
        &dir.join("concentrations.csv"),
        &table_csv(&labels, res.t(), &rows),
    )?;
    let (labels, rows) = res.reaction_rates();
    write(
        &dir.join("reaction_rates.csv"),
        &table_csv(&labels, res.t(), &rows),
    )?;
    let integrals = res.integrate_species(0.0, f64::INFINITY);
    write(
        &dir.join("species_integrals.csv"),
        &series_csv("species", &integrals),
    )?;
    for sp in production {
        let v = res.integrate_production_of(sp, None, None, true, true)?;
        write(
            &dir.join(format!("production_{sp}.csv")),
            &series_csv("reaction", &v),
        )?;
    }
    for sp in consumption {
        let v = res.integrate_consumption_of(sp, None, None, true)?;
        write(
            &dir.join(format!("consumption_{sp}.csv")),
            &series_csv("reaction", &v),
        )?;
    }
    let json = serde_json::to_string_pretty(meta).map_err(|e| RadkinError::Value(e.to_string()))?;
    write(&dir.join("run.json"), &json)
}
