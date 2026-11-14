//! Files of a voxel run.
//!
//! | file | content |
//! |---|---|
//! | `yields.csv` | `t`, `N_<sp>` (molecules in the domain) and, when seeds give `energy_ev`, `G_<sp>` = (N - N(t_start)) / (E / 100 eV) with E the energy of the seeds injected so far [#/100 eV] |
//! | `snapshots/t<k>_<sp>.npy` | concentration field [mol/L], float64, C order, shape `(nz, ny, nx)` |
//! | `snapshots.json` | grid, species, snapshot times and file names, run statistics |

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use serde::Serialize;

use crate::consts::AVOGADRO;
use crate::error::{RadkinError, Result};
use crate::voxel::grid::VoxelGrid;
use crate::voxel::sim::VoxelOutput;

/// `.npy` (format 1.0) bytes of a little-endian f64 array.
pub fn npy_bytes(data: &[f64], shape: &[usize]) -> Vec<u8> {
    let dims: Vec<String> = shape.iter().map(|d| d.to_string()).collect();
    let shape_str = if dims.len() == 1 {
        format!("({},)", dims[0])
    } else {
        format!("({})", dims.join(", "))
    };
    let mut header = format!("{{'descr': '<f8', 'fortran_order': False, 'shape': {shape_str}, }}");
    // magic(6) + version(2) + len(2) + header + '\n' padded to 64 bytes.
    let total = 10 + header.len() + 1;
    header.push_str(&" ".repeat((64 - total % 64) % 64));
    header.push('\n');
    let mut out = Vec::with_capacity(10 + header.len() + 8 * data.len());
    out.extend_from_slice(b"\x93NUMPY\x01\x00");
    out.extend_from_slice(&(header.len() as u16).to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    for v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

#[derive(Serialize)]
struct SnapshotIndex<'a> {
    dims: [usize; 3],
    spacing_nm: [f64; 3],
    origin_nm: [f64; 3],
    species: &'a [String],
    unit: &'static str,
    layout: &'static str,
    snapshots: Vec<SnapshotEntry>,
    n_steps: usize,
    n_rejected: usize,
}

#[derive(Serialize)]
struct SnapshotEntry {
    t: f64,
    files: Vec<String>,
}

fn write(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(|e| RadkinError::io(path, e))?;
    f.write_all(bytes).map_err(|e| RadkinError::io(path, e))
}

/// Molecules per µmol/L summed over voxels.
fn molecules_per_unit(grid: &VoxelGrid) -> f64 {
    1e-6 * AVOGADRO * grid.voxel_volume_l()
}

/// `yields.csv` text.
pub fn yields_csv(out: &VoxelOutput, grid: &VoxelGrid) -> String {
    let k = molecules_per_unit(grid);
    let with_g = !out.seed_energies.is_empty();
    let mut s = String::from("t");
    for sp in &out.species {
        let _ = write!(s, ",N_{sp}");
    }
    if with_g {
        for sp in &out.species {
            let _ = write!(s, ",G_{sp}");
        }
    }
    s.push('\n');
    let n0 = out
        .totals
        .first()
        .map(|(_, v)| v.clone())
        .unwrap_or_default();
    for (i, (t, sums)) in out.totals.iter().enumerate() {
        // A step landing on an injection time is recorded before and after the
        // injection: keep only the state after it.
        if out.totals.get(i + 1).is_some_and(|(next, _)| next == t) {
            continue;
        }
        let _ = write!(s, "{t:e}");
        for v in sums {
            let _ = write!(s, ",{:e}", v * k);
        }
        if with_g {
            // Net of the background chemistry when it is tracked (the
            // background may drift, e.g. a pH not exactly at equilibrium).
            let base = out
                .background_totals
                .get(i)
                .filter(|_| out.background_totals.len() == out.totals.len());
            let base = base.unwrap_or(&n0);
            // Seeds injected at `t` are already in `sums` (recorded after the injection).
            let e: f64 = out
                .seed_energies
                .iter()
                .filter(|(ts, _)| ts <= t)
                .map(|(_, e)| e)
                .sum();
            for (v, v0) in sums.iter().zip(base) {
                if e > 0.0 {
                    let _ = write!(s, ",{:e}", (v - v0) * k / (e / 100.0));
                } else {
                    s.push_str(",nan");
                }
            }
        }
        s.push('\n');
    }
    s
}

/// Write `yields.csv`, `snapshots/` and `snapshots.json` into `dir`.
pub fn write_voxel_dir(dir: &Path, out: &VoxelOutput, grid: &VoxelGrid) -> Result<()> {
    let snap_dir = dir.join("snapshots");
    std::fs::create_dir_all(&snap_dir).map_err(|e| RadkinError::io(&snap_dir, e))?;
    write(&dir.join("yields.csv"), yields_csv(out, grid).as_bytes())?;
    let [nx, ny, nz] = grid.dims;
    let mut entries = Vec::new();
    for (k, snap) in out.snapshots.iter().enumerate() {
        let mut files = Vec::new();
        for (i, sp) in out.species.iter().enumerate() {
            let name = format!("t{k:03}_{sp}.npy");
            let molar: Vec<f64> = snap.state.field(i).iter().map(|v| v * 1e-6).collect();
            write(&snap_dir.join(&name), &npy_bytes(&molar, &[nz, ny, nx]))?;
            files.push(format!("snapshots/{name}"));
        }
        entries.push(SnapshotEntry { t: snap.t, files });
    }
    let index = SnapshotIndex {
        dims: grid.dims,
        spacing_nm: grid.spacing.map(|x| x * 1e9),
        origin_nm: grid.origin.map(|x| x * 1e9),
        species: &out.species,
        unit: "mol/L",
        layout: "C order, shape (nz, ny, nx)",
        snapshots: entries,
        n_steps: out.n_steps,
        n_rejected: out.n_rejected,
    };
    let json = serde_json::to_string_pretty(&index)
        .map_err(|e| RadkinError::Value(format!("snapshots.json: {e}")))?;
    write(&dir.join("snapshots.json"), json.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn yields_g_values_use_injected_energy() {
        let grid = VoxelGrid::new([1, 1, 1], [10e-9; 3], [0.0; 3]).unwrap();
        let k = molecules_per_unit(&grid);
        let out = VoxelOutput {
            species: vec!["A".into()],
            // 0 molecules, then 10 after a 100 eV seed, then 5.
            totals: vec![
                (0.0, vec![0.0]),
                (1e-12, vec![10.0 / k]),
                (1e-9, vec![5.0 / k]),
            ],
            background_totals: Vec::new(),
            snapshots: Vec::new(),
            seed_energies: vec![(1e-12, 100.0)],
            n_steps: 2,
            n_rejected: 0,
            level_history: Vec::new(),
        };
        let csv = yields_csv(&out, &grid);
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 4, "{csv}");
        assert_eq!(lines[0], "t,N_A,G_A");
        assert!(lines[1].ends_with(",nan"));
        let g = |l: &str| l.rsplit(',').next().unwrap().parse::<f64>().unwrap();
        assert!((g(lines[2]) - 10.0).abs() < 1e-9);
        assert!((g(lines[3]) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn yields_g_values_are_net_of_the_background() {
        let grid = VoxelGrid::new([1, 1, 1], [10e-9; 3], [0.0; 3]).unwrap();
        let k = molecules_per_unit(&grid);
        // The background alone goes 100 -> 100 -> 104 molecules; the seed adds 10.
        let out = VoxelOutput {
            species: vec!["A".into()],
            totals: vec![
                (0.0, vec![100.0 / k]),
                (1e-12, vec![110.0 / k]),
                (1e-9, vec![112.0 / k]),
            ],
            background_totals: vec![vec![100.0 / k], vec![100.0 / k], vec![104.0 / k]],
            snapshots: Vec::new(),
            seed_energies: vec![(1e-12, 100.0)],
            n_steps: 2,
            n_rejected: 0,
            level_history: Vec::new(),
        };
        let csv = yields_csv(&out, &grid);
        let g = |l: &str| l.rsplit(',').next().unwrap().parse::<f64>().unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert!((g(lines[2]) - 10.0).abs() < 1e-9, "{csv}");
        assert!((g(lines[3]) - 8.0).abs() < 1e-9, "{csv}");
    }

    #[test]
    fn npy_header_is_aligned_and_parsable() {
        let b = npy_bytes(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[1, 2, 3]);
        assert_eq!(&b[..8], b"\x93NUMPY\x01\x00");
        let hlen = u16::from_le_bytes([b[8], b[9]]) as usize;
        assert_eq!((10 + hlen) % 64, 0);
        let header = std::str::from_utf8(&b[10..10 + hlen]).unwrap();
        assert!(header.contains("'shape': (1, 2, 3)"), "{header}");
        assert!(header.ends_with('\n'));
        assert_eq!(b.len(), 10 + hlen + 48);
        assert_eq!(
            f64::from_le_bytes(b[10 + hlen..18 + hlen].try_into().unwrap()),
            1.0
        );
    }
}
