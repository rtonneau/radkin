//! Files of an IRT run.
//!
//! | file | content |
//! |---|---|
//! | `yields.csv` | `t`, `N_<sp>` (mean particles), `G_<sp>` when seeds give `energy_ev` (`N / (E / 100 eV)`, E the energy of the seeds injected so far), then `N_<sp>_std` and `G_<sp>_std` (standard deviation over the trajectories): the columns of the RDME `yields.csv` |
//! | `irt.json` | settings, event counts, firings of each reaction, molecules of the background consumed, particles that left the box, warnings |

use std::fmt::Write as _;
use std::path::Path;

use indexmap::IndexMap;
use serde::Serialize;

use crate::error::{RadkinError, Result};
use crate::irt::sim::{IrtOutput, IrtSimulation};

fn header(s: &mut String, species: &[String], prefix: &str, suffix: &str) {
    for sp in species {
        let _ = write!(s, ",{prefix}_{sp}{suffix}");
    }
}

/// `values` as CSV fields, each divided by `divisor` (`nan` if it is not > 0).
fn fields(s: &mut String, values: &[f64], divisor: f64) {
    for v in values {
        if divisor > 0.0 {
            let _ = write!(s, ",{:e}", v / divisor);
        } else {
            s.push_str(",nan");
        }
    }
}

/// `yields.csv` text.
pub fn irt_yields_csv(out: &IrtOutput) -> String {
    let with_g = !out.seed_energies.is_empty();
    let mut s = String::from("t");
    header(&mut s, &out.species, "N", "");
    if with_g {
        header(&mut s, &out.species, "G", "");
    }
    header(&mut s, &out.species, "N", "_std");
    if with_g {
        header(&mut s, &out.species, "G", "_std");
    }
    s.push('\n');
    for (k, t) in out.record_times.iter().enumerate() {
        let _ = write!(s, "{t:e}");
        // Seeds injected at `t` are in the totals of `t`. G is per 100 eV.
        let energy: f64 = out
            .seed_energies
            .iter()
            .filter(|(ts, _)| ts <= t)
            .map(|(_, e)| e)
            .sum();
        fields(&mut s, &out.totals_mean[k], 1.0);
        if with_g {
            fields(&mut s, &out.totals_mean[k], energy / 100.0);
        }
        fields(&mut s, &out.totals_std[k], 1.0);
        if with_g {
            fields(&mut s, &out.totals_std[k], energy / 100.0);
        }
        s.push('\n');
    }
    s
}

#[derive(Serialize)]
struct Meta<'a> {
    solver: &'static str,
    resampling: &'static str,
    t_start: f64,
    t_end: f64,
    trajectories: usize,
    seed: u64,
    start_contact: bool,
    boundary: &'static str,
    box_nm: Option<[[f64; 3]; 2]>,
    cutoff_probability: f64,
    species: &'a [String],
    /// Molecules injected by the seeds, per trajectory.
    seed_molecules: usize,
    seed_energy_ev: f64,
    /// Event counts summed over the trajectories.
    pair_reactions: u64,
    contact_reactions: u64,
    /// Reactions of one particle (first order, with a fixed or a background
    /// species).
    molecule_reactions: u64,
    background_reactions: bool,
    /// Molecules of each background species these reactions consumed,
    /// summed over the trajectories. The background is not debited.
    background_consumed: IndexMap<&'a str, u64>,
    /// Concentration of each background species [mol/L].
    background: IndexMap<&'a str, f64>,
    sampling_passes: u64,
    /// Moves made to look again at a pair inside its effective radius.
    look_again_moves: u64,
    pairs_sampled: u64,
    /// Firings of each reaction that fired, summed over the trajectories.
    reactions: IndexMap<&'a str, u64>,
    /// Particles that left the box, summed over the trajectories.
    absorbed: IndexMap<&'a str, u64>,
    pair_reaction_types: Vec<PairType<'a>>,
    warnings: &'a [String],
    exec_time: f64,
    radkin_version: &'static str,
}

#[derive(Serialize)]
struct PairType<'a> {
    reaction: &'a str,
    irt_type: &'static str,
    typed_by: &'static str,
    r_eff_nm: f64,
}

/// `irt.json` text.
pub fn irt_meta_json(sim: &IrtSimulation, out: &IrtOutput) -> Result<String> {
    use crate::irt::config::BoundaryKind;
    use crate::irt::model::TypeSource;
    let meta = Meta {
        solver: "irt",
        resampling: "full",
        t_start: sim.t_start,
        t_end: sim.config.t_end,
        trajectories: out.trajectories,
        seed: out.seed,
        start_contact: sim.config.start_contact,
        boundary: match sim.config.boundary {
            BoundaryKind::None => "none",
            BoundaryKind::Outflow => "outflow",
        },
        box_nm: sim
            .config
            .domain
            .map(|(lo, hi)| [lo.map(|x| x * 1e9), hi.map(|x| x * 1e9)]),
        cutoff_probability: sim.config.cutoff_probability,
        species: &out.species,
        seed_molecules: sim.seeds.iter().map(|s| s.molecules.len()).sum(),
        seed_energy_ev: out.seed_energies.iter().map(|(_, e)| e).sum(),
        pair_reactions: out.stats.pair_reactions,
        contact_reactions: out.stats.contact_reactions,
        molecule_reactions: out.stats.molecule_reactions,
        background_reactions: sim.config.background_reactions,
        background_consumed: out
            .species
            .iter()
            .zip(&out.stats.background_consumed)
            .filter(|(_, n)| **n > 0)
            .map(|(s, n)| (s.as_str(), *n))
            .collect(),
        background: out
            .species
            .iter()
            .zip(&sim.model.background)
            .filter(|(_, c)| **c > 0.0)
            .map(|(s, c)| (s.as_str(), *c))
            .collect(),
        sampling_passes: out.stats.passes,
        look_again_moves: out.stats.moves,
        pairs_sampled: out.stats.pairs_sampled,
        reactions: out
            .reactions
            .iter()
            .zip(&out.stats.by_reaction)
            .filter(|(_, n)| **n > 0)
            .map(|(r, n)| (r.as_str(), *n))
            .collect(),
        absorbed: out
            .species
            .iter()
            .zip(&out.stats.absorbed)
            .filter(|(_, n)| **n > 0)
            .map(|(s, n)| (s.as_str(), *n))
            .collect(),
        pair_reaction_types: sim
            .model
            .pairs
            .iter()
            .map(|p| PairType {
                reaction: &out.reactions[p.index],
                irt_type: p.kind.name(),
                typed_by: match p.source {
                    TypeSource::Network => "network",
                    TypeSource::Rule => "rule",
                    TypeSource::Fallback => "fallback",
                },
                r_eff_nm: p.r_eff * 1e9,
            })
            .collect(),
        warnings: &out.warnings,
        exec_time: out.exec_time,
        radkin_version: env!("CARGO_PKG_VERSION"),
    };
    serde_json::to_string_pretty(&meta).map_err(|e| RadkinError::Value(e.to_string()))
}

/// Write `yields.csv` and `irt.json` into `dir`.
pub fn write_irt_dir(dir: &Path, sim: &IrtSimulation, out: &IrtOutput) -> Result<()> {
    std::fs::create_dir_all(dir).map_err(|e| RadkinError::io(dir, e))?;
    let write = |name: &str, text: String| {
        let path = dir.join(name);
        std::fs::write(&path, text).map_err(|e| RadkinError::io(&path, e))
    };
    write("yields.csv", irt_yields_csv(out))?;
    write("irt.json", irt_meta_json(sim, out)?)
}
