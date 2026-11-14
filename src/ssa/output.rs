//! Output files of a unit-cell SSA run.
//!
//! | file | content |
//! |---|---|
//! | `concentrations.csv` | `t`, ensemble mean of the dynamic species, the fixed species, then `<species>_std` of the dynamic species (µmol/L) |
//! | `run.json` | inputs, G-values and the SSA settings |
//!
//! The first columns are those of the LSODA `concentrations.csv`.

use std::path::Path;

use indexmap::IndexMap;
use serde::Serialize;

use crate::error::{RadkinError, Result};
use crate::output::table_csv;
use crate::ssa::ensemble::SsaResult;

/// Run metadata written to `run.json`.
#[derive(Debug, Clone, Serialize)]
pub struct SsaRunMeta {
    pub config: String,
    pub solver: String,
    pub t_span: (f64, f64),
    pub beam_label: Option<String>,
    pub beam_type: Option<String>,
    pub beam_params: IndexMap<String, f64>,
    #[serde(rename = "LET")]
    pub let_: Option<f64>,
    pub temperature: f64,
    pub g_values: IndexMap<String, f64>,
    /// Volume [L].
    pub volume: f64,
    /// One molecule in the volume [µmol/L].
    pub molecule_concentration: f64,
    pub trajectories: usize,
    pub seed: u64,
    pub tau_leap: bool,
    pub n_points: usize,
    /// Events fired, summed over the trajectories.
    pub n_events: u64,
    pub exec_time: f64,
    pub radkin_version: String,
    pub species_labels: Vec<String>,
    /// `as_label() -> str(reaction)`.
    pub reaction_strings: IndexMap<String, String>,
}

impl SsaRunMeta {
    pub fn from_result(
        res: &SsaResult,
        config: &str,
        t_span: (f64, f64),
        beam_label: Option<&str>,
    ) -> Self {
        let beam_params = res
            .beam
            .as_ref()
            .map(|b| {
                b.as_table()
                    .into_iter()
                    .filter_map(|(k, v)| v.as_float().map(|f| (k, f)))
                    .collect()
            })
            .unwrap_or_default();
        let mut species_labels = res.env.species.sorted_dyn_labels();
        species_labels.extend(res.env.species.sorted_cst_labels());
        Self {
            config: config.to_string(),
            solver: "ssa".into(),
            t_span,
            beam_label: res.beam.as_ref().and(beam_label.map(String::from)),
            beam_type: res.beam.as_ref().map(|b| b.beam_type().to_string()),
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
            volume: res.volume,
            molecule_concentration: res.molecule(),
            trajectories: res.trajectories,
            seed: res.seed,
            tau_leap: res.tau_leap,
            n_points: res.t.len(),
            n_events: res.n_events,
            exec_time: res.exec_time,
            radkin_version: env!("CARGO_PKG_VERSION").to_string(),
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

/// `concentrations.csv` text: the LSODA columns (mean), then the `_std` columns.
pub fn concentrations_csv(res: &SsaResult) -> String {
    let (mut labels, mut rows) = res.concentrations();
    labels.extend(
        res.env
            .species
            .sorted_dyn_labels()
            .into_iter()
            .map(|l| format!("{l}_std")),
    );
    for (row, std) in rows.iter_mut().zip(&res.std) {
        row.extend(std);
    }
    table_csv(&labels, &res.t, &rows)
}

/// Write `concentrations.csv` and `run.json` into `dir`.
pub fn write_ssa_dir(dir: &Path, res: &SsaResult, meta: &SsaRunMeta) -> Result<()> {
    std::fs::create_dir_all(dir).map_err(|e| RadkinError::io(dir, e))?;
    let write = |name: &str, text: &str| {
        let path = dir.join(name);
        std::fs::write(&path, text).map_err(|e| RadkinError::io(&path, e))
    };
    write("concentrations.csv", &concentrations_csv(res))?;
    let json = serde_json::to_string_pretty(meta).map_err(|e| RadkinError::Value(e.to_string()))?;
    write("run.json", &json)
}
