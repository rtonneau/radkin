//! Files of an RDME voxel run: the layout of the deterministic voxel run
//! (ensemble mean), plus the standard deviation over the trajectories.
//!
//! | file | content |
//! |---|---|
//! | `yields.csv` | `t`, `N_<sp>` (mean molecules in the domain), `G_<sp>` when seeds give `energy_ev` (mean of `(N - N(t_start)) / (E / 100 eV)`, E the energy of the seeds injected so far), then `N_<sp>_std` and `G_<sp>_std` |
//! | | With `[ssa] background_run` (default) the `G` and `G_std` are net of the seedless background run's totals (means subtracted, standard deviations added in quadrature) instead of the totals at `t_start`; `N` stays gross |
//! | `snapshots/t<k>_<sp>.npy` | mean concentration field [mol/L], float64, C order, shape `(nz, ny, nx)` |
//! | `snapshots/t<k>_<sp>_std.npy` | standard deviation of the field over the trajectories [mol/L] |
//! | `snapshots.json` | grid, species, snapshot times and file names, ensemble settings and event counts |

use std::fmt::Write as _;
use std::path::Path;

use serde::Serialize;

use crate::error::{RadkinError, Result};
use crate::ssa::voxel_sim::RdmeOutput;
use crate::voxel::grid::VoxelGrid;
use crate::voxel::output::npy_bytes;

#[derive(Serialize)]
struct SnapshotIndex<'a> {
    dims: [usize; 3],
    spacing_nm: [f64; 3],
    origin_nm: [f64; 3],
    species: &'a [String],
    unit: &'static str,
    layout: &'static str,
    snapshots: Vec<SnapshotEntry>,
    solver: &'static str,
    trajectories: usize,
    seed: u64,
    tau_leap: bool,
    rate_correction: &'static str,
    /// The seedless run subtracted from `yields.csv`, with its noise.
    background_run: Option<BackgroundIndex<'a>>,
    /// Events summed over the trajectories.
    n_reactions: u64,
    n_hops: u64,
    /// Species level of every species, in `species` order.
    species_levels: &'a [u32],
    exec_time: f64,
}

#[derive(Serialize)]
struct BackgroundIndex<'a> {
    /// Seed of the background run, derived from `seed`.
    seed: u64,
    /// Standard error of the mean change of every species (`species` order)
    /// in the background run over the whole run [molecules]: the noise the
    /// subtraction adds to `yields.csv`.
    noise: &'a [f64],
    n_reactions: u64,
    n_hops: u64,
}

#[derive(Serialize)]
struct SnapshotEntry {
    t: f64,
    /// Ensemble mean, one file per species.
    files: Vec<String>,
    /// Standard deviation over the trajectories, one file per species.
    std_files: Vec<String>,
}

fn write(path: &Path, bytes: &[u8]) -> Result<()> {
    std::fs::write(path, bytes).map_err(|e| RadkinError::io(path, e))
}

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
pub fn rdme_yields_csv(out: &RdmeOutput) -> String {
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
            fields(&mut s, &out.net_mean[k], energy / 100.0);
        }
        fields(&mut s, &out.totals_std[k], 1.0);
        if with_g {
            fields(&mut s, &out.net_std[k], energy / 100.0);
        }
        s.push('\n');
    }
    s
}

/// Write `yields.csv`, `snapshots/` and `snapshots.json` into `dir`.
pub fn write_rdme_dir(dir: &Path, out: &RdmeOutput, grid: &VoxelGrid) -> Result<()> {
    let snap_dir = dir.join("snapshots");
    std::fs::create_dir_all(&snap_dir).map_err(|e| RadkinError::io(&snap_dir, e))?;
    write(&dir.join("yields.csv"), rdme_yields_csv(out).as_bytes())?;
    let [nx, ny, nz] = grid.dims;
    let mut entries = Vec::new();
    for (k, t) in out.snapshot_times.iter().enumerate() {
        let (mut files, mut std_files) = (Vec::new(), Vec::new());
        for (i, sp) in out.species.iter().enumerate() {
            for (state, suffix, list) in [
                (&out.mean[k], "", &mut files),
                (&out.std[k], "_std", &mut std_files),
            ] {
                let name = format!("t{k:03}_{sp}{suffix}.npy");
                let molar: Vec<f64> = state.field(i).iter().map(|v| v * 1e-6).collect();
                write(&snap_dir.join(&name), &npy_bytes(&molar, &[nz, ny, nx]))?;
                list.push(format!("snapshots/{name}"));
            }
        }
        entries.push(SnapshotEntry {
            t: *t,
            files,
            std_files,
        });
    }
    let index = SnapshotIndex {
        dims: grid.dims,
        spacing_nm: grid.spacing.map(|x| x * 1e9),
        origin_nm: grid.origin.map(|x| x * 1e9),
        species: &out.species,
        unit: "mol/L",
        layout: "C order, shape (nz, ny, nx)",
        snapshots: entries,
        solver: "ssa",
        trajectories: out.trajectories,
        seed: out.seed,
        tau_leap: out.tau_leap,
        rate_correction: out.rate_correction,
        background_run: out.background.as_ref().map(|b| BackgroundIndex {
            seed: b.seed,
            noise: &b.noise,
            n_reactions: b.n_reactions,
            n_hops: b.n_hops,
        }),
        n_reactions: out.n_reactions,
        n_hops: out.n_hops,
        species_levels: &out.species_levels,
        exec_time: out.exec_time,
    };
    let json = serde_json::to_string_pretty(&index)
        .map_err(|e| RadkinError::Value(format!("snapshots.json: {e}")))?;
    write(&dir.join("snapshots.json"), json.as_bytes())
}
