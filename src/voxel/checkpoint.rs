//! Voxel run checkpoints: the state of a [`VoxelRun`](crate::voxel::VoxelRun)
//! at one time, saved as a NumPy `.npz` archive (uncompressed zip of `.npy`
//! arrays plus a JSON entry), so `np.load` reads it directly.
//!
//! | Entry | Content |
//! |---|---|
//! | `state.npy` | `<f8`, shape `(n_species, nz, ny, nx)` [µmol/L] |
//! | `background.npy` | `<f8`, shape `(n_species,)`, only with a background |
//! | `time.npy` | `<f8` `[t, dt]` [s] |
//! | `grid.npy` | `<f8` `[dx, dy, dz, x0, y0, z0]` [m] |
//! | `seeds.npy` | `<f8` `[time_0, ..]` of the injected seeds [s] |
//! | `meta.json` | version, dims, species, seed fingerprints, config summary |
//!
//! Floats are stored in `.npy` arrays so they round-trip bit for bit.
//! Writing goes to a temporary file renamed over the target, so a crash
//! never leaves a half-written checkpoint.

use std::io::{Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{RadkinError, Result};
use crate::voxel::grid::{VoxelGrid, VoxelState};

/// Checkpoint format version written in `meta.json`.
pub const CHECKPOINT_VERSION: u32 = 1;

/// A seed already injected when the checkpoint was taken.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeedRecord {
    /// Injection time [s].
    pub time: f64,
    /// Hash of the injected amounts (see [`SeedRecord::fingerprint_of`]).
    pub fingerprint: u64,
}

impl SeedRecord {
    /// FNV-1a hash of the bit patterns of `amount`.
    pub fn fingerprint_of(amount: &VoxelState) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for v in &amount.data {
            for b in v.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        h
    }
}

/// The state of a voxel run at time `t`.
#[derive(Debug, Clone, PartialEq)]
pub struct Checkpoint {
    pub t: f64,
    /// Splitting step size to resume with [s].
    pub dt: f64,
    pub grid: VoxelGrid,
    /// Dynamic species labels in state order.
    pub species: Vec<String>,
    pub state: VoxelState,
    /// Background vector [µmol/L], when the run tracks one.
    pub background: Option<Vec<f64>>,
    pub seeds: Vec<SeedRecord>,
    /// Free-form summary of the run configuration, compared on resume.
    pub config: serde_json::Value,
}

#[derive(Serialize, Deserialize)]
struct Meta {
    version: u32,
    dims: [usize; 3],
    species: Vec<String>,
    seed_fingerprints: Vec<u64>,
    config: serde_json::Value,
}

fn ckpt_err(path: &Path, msg: impl std::fmt::Display) -> RadkinError {
    RadkinError::Checkpoint(format!("{}: {msg}", path.display()))
}

/// `.npy` (format 1.0) bytes of a little-endian f64 array.
fn npy_bytes(shape: &[usize], data: &[f64]) -> Vec<u8> {
    debug_assert_eq!(shape.iter().product::<usize>(), data.len());
    let dims: Vec<String> = shape.iter().map(|d| d.to_string()).collect();
    let shape_txt = if dims.len() == 1 {
        format!("({},)", dims[0])
    } else {
        format!("({})", dims.join(", "))
    };
    let mut header = format!("{{'descr': '<f8', 'fortran_order': False, 'shape': {shape_txt}, }}");
    // Magic (6) + version (2) + length (2) + header + '\n', padded to 64.
    let unpadded = 10 + header.len() + 1;
    header.push_str(&" ".repeat((64 - unpadded % 64) % 64));
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

/// Parse `.npy` bytes written by [`npy_bytes`] (or numpy, `<f8`, C order).
fn parse_npy(bytes: &[u8]) -> std::result::Result<(Vec<usize>, Vec<f64>), String> {
    if bytes.len() < 10 || &bytes[..6] != b"\x93NUMPY" {
        return Err("not a .npy array".into());
    }
    let (hlen, start) = match bytes[6] {
        1 => (u16::from_le_bytes([bytes[8], bytes[9]]) as usize, 10),
        2 | 3 if bytes.len() >= 12 => (
            u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize,
            12,
        ),
        v => return Err(format!("unsupported .npy version {v}")),
    };
    let header = bytes
        .get(start..start + hlen)
        .and_then(|h| std::str::from_utf8(h).ok())
        .ok_or("truncated .npy header")?;
    if !header.contains("'descr': '<f8'") {
        return Err(format!("expected little-endian f64 array, header {header}"));
    }
    if header.contains("'fortran_order': True") {
        return Err("fortran-ordered arrays are not supported".into());
    }
    let shape_txt = header
        .split("'shape':")
        .nth(1)
        .and_then(|r| r.split('(').nth(1))
        .and_then(|r| r.split(')').next())
        .ok_or_else(|| format!("no shape in header {header}"))?;
    let shape = shape_txt
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .map_err(|e| format!("bad shape '{s}': {e}"))
        })
        .collect::<std::result::Result<Vec<usize>, String>>()?;
    let n: usize = shape.iter().product();
    let body = &bytes[start + hlen..];
    if body.len() != 8 * n {
        return Err(format!("expected {n} values, found {} bytes", body.len()));
    }
    let data = body
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| f64::from_le_bytes(*c))
        .collect();
    Ok((shape, data))
}

impl Checkpoint {
    /// Write the checkpoint to `path` (temporary file, then rename).
    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let [nx, ny, nz] = self.grid.dims;
        if self.state.n_species != self.species.len() || self.state.n_voxels != nx * ny * nz {
            return Err(ckpt_err(path, "state does not match the grid and species"));
        }
        let meta = Meta {
            version: CHECKPOINT_VERSION,
            dims: self.grid.dims,
            species: self.species.clone(),
            seed_fingerprints: self.seeds.iter().map(|s| s.fingerprint).collect(),
            config: self.config.clone(),
        };
        let meta = serde_json::to_vec_pretty(&meta).map_err(|e| ckpt_err(path, e))?;
        let g = &self.grid;
        let seed_times: Vec<f64> = self.seeds.iter().map(|s| s.time).collect();
        let mut entries: Vec<(&str, Vec<u8>)> = vec![
            (
                "state.npy",
                npy_bytes(&[self.state.n_species, nz, ny, nx], &self.state.data),
            ),
            ("time.npy", npy_bytes(&[2], &[self.t, self.dt])),
            (
                "grid.npy",
                npy_bytes(
                    &[6],
                    &[
                        g.spacing[0],
                        g.spacing[1],
                        g.spacing[2],
                        g.origin[0],
                        g.origin[1],
                        g.origin[2],
                    ],
                ),
            ),
            ("seeds.npy", npy_bytes(&[seed_times.len()], &seed_times)),
            ("meta.json", meta),
        ];
        if let Some(b) = &self.background {
            entries.push(("background.npy", npy_bytes(&[b.len()], b)));
        }

        let mut tmp = path.as_os_str().to_owned();
        tmp.push(".tmp");
        let tmp = std::path::PathBuf::from(tmp);
        let file = std::fs::File::create(&tmp).map_err(|e| RadkinError::io(&tmp, e))?;
        let mut zip = zip::ZipWriter::new(std::io::BufWriter::new(file));
        let opts = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .large_file(entries.iter().any(|(_, b)| b.len() >= u32::MAX as usize));
        for (name, bytes) in &entries {
            zip.start_file(*name, opts).map_err(|e| ckpt_err(&tmp, e))?;
            zip.write_all(bytes).map_err(|e| RadkinError::io(&tmp, e))?;
        }
        let mut w = zip.finish().map_err(|e| ckpt_err(&tmp, e))?;
        w.flush().map_err(|e| RadkinError::io(&tmp, e))?;
        w.get_ref()
            .sync_all()
            .map_err(|e| RadkinError::io(&tmp, e))?;
        drop(w);
        std::fs::rename(&tmp, path).map_err(|e| RadkinError::io(path, e))
    }

    /// Read a checkpoint written by [`Checkpoint::write`].
    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = std::fs::File::open(path).map_err(|e| RadkinError::io(path, e))?;
        let mut zip =
            zip::ZipArchive::new(std::io::BufReader::new(file)).map_err(|e| ckpt_err(path, e))?;
        let mut entry = |name: &str| -> Result<Option<Vec<u8>>> {
            let mut f = match zip.by_name(name) {
                Ok(f) => f,
                Err(zip::result::ZipError::FileNotFound) => return Ok(None),
                Err(e) => return Err(ckpt_err(path, e)),
            };
            let mut buf = Vec::with_capacity(f.size() as usize);
            f.read_to_end(&mut buf)
                .map_err(|e| RadkinError::io(path, e))?;
            Ok(Some(buf))
        };
        let required = |v: Option<Vec<u8>>, name: &str| {
            v.ok_or_else(|| ckpt_err(path, format!("missing entry {name}")))
        };
        let array = |bytes: Vec<u8>, name: &str| {
            parse_npy(&bytes).map_err(|e| ckpt_err(path, format!("{name}: {e}")))
        };

        let meta: Meta = serde_json::from_slice(&required(entry("meta.json")?, "meta.json")?)
            .map_err(|e| ckpt_err(path, format!("meta.json: {e}")))?;
        if meta.version != CHECKPOINT_VERSION {
            return Err(ckpt_err(
                path,
                format!(
                    "format version {} (this radkin reads {CHECKPOINT_VERSION})",
                    meta.version
                ),
            ));
        }
        let (_, time) = array(required(entry("time.npy")?, "time.npy")?, "time.npy")?;
        let (_, g) = array(required(entry("grid.npy")?, "grid.npy")?, "grid.npy")?;
        let (_, seed_times) = array(required(entry("seeds.npy")?, "seeds.npy")?, "seeds.npy")?;
        let (shape, data) = array(required(entry("state.npy")?, "state.npy")?, "state.npy")?;
        let background = entry("background.npy")?
            .map(|b| array(b, "background.npy").map(|(_, v)| v))
            .transpose()?;

        if time.len() != 2 || g.len() != 6 || seed_times.len() != meta.seed_fingerprints.len() {
            return Err(ckpt_err(path, "inconsistent time, grid or seed entries"));
        }
        let [nx, ny, nz] = meta.dims;
        let n_species = meta.species.len();
        if shape != [n_species, nz, ny, nx] {
            return Err(ckpt_err(
                path,
                format!(
                    "state shape {shape:?} does not match {n_species} species on {:?}",
                    meta.dims
                ),
            ));
        }
        if background.as_ref().is_some_and(|b| b.len() != n_species) {
            return Err(ckpt_err(path, "background length differs from species"));
        }
        let grid = VoxelGrid::new(meta.dims, [g[0], g[1], g[2]], [g[3], g[4], g[5]])
            .map_err(|e| ckpt_err(path, e))?;
        Ok(Self {
            t: time[0],
            dt: time[1],
            state: VoxelState {
                n_species,
                n_voxels: nx * ny * nz,
                data,
            },
            grid,
            species: meta.species,
            background,
            seeds: seed_times
                .into_iter()
                .zip(meta.seed_fingerprints)
                .map(|(time, fingerprint)| SeedRecord { time, fingerprint })
                .collect(),
            config: meta.config,
        })
    }

    /// Error unless the checkpoint has the same grid (dims, spacing, origin)
    /// and the same species in the same order.
    pub fn check_compatible(&self, grid: &VoxelGrid, species: &[String]) -> Result<()> {
        if self.grid != *grid {
            return Err(RadkinError::Checkpoint(format!(
                "grid differs: checkpoint {:?} spacing {:?} origin {:?}, config {:?} spacing {:?} origin {:?}",
                self.grid.dims,
                self.grid.spacing,
                self.grid.origin,
                grid.dims,
                grid.spacing,
                grid.origin
            )));
        }
        if self.species != species {
            return Err(RadkinError::Checkpoint(format!(
                "species differ: checkpoint {:?}, config {species:?}",
                self.species
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(background: bool) -> Checkpoint {
        let grid = VoxelGrid::new([3, 2, 2], [1e-9, 2e-9, 3e-9], [-1.5e-9, 0.0, 0.1]).unwrap();
        let species = vec!["e_aq".to_string(), "°OH".to_string()];
        let mut state = VoxelState::zeros(2, grid.n_voxels());
        for (i, v) in state.data.iter_mut().enumerate() {
            *v = 0.1 + i as f64 / 3.0;
        }
        Checkpoint {
            t: 1.234_567_890_123e-9,
            dt: 3.3e-12,
            grid,
            species,
            state,
            background: background.then(|| vec![1e-7, 0.1 / 3.0]),
            seeds: vec![SeedRecord {
                time: 1e-12,
                fingerprint: u64::MAX - 7,
            }],
            config: serde_json::json!({"t_end": 1e-6}),
        }
    }

    fn tmp_path(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("radkin-checkpoint-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    #[test]
    fn checkpoint_round_trip_is_exact() {
        for bg in [false, true] {
            let c = sample(bg);
            let p = tmp_path(&format!("round_trip_{bg}.npz"));
            c.write(&p).unwrap();
            assert_eq!(Checkpoint::read(&p).unwrap(), c);
            let mut tmp = p.as_os_str().to_owned();
            tmp.push(".tmp");
            assert!(!std::path::Path::new(&tmp).exists());
        }
    }

    #[test]
    fn checkpoint_npy_header_is_numpy_format() {
        let b = npy_bytes(&[2, 3], &[0.0; 6]);
        assert_eq!(&b[..8], b"\x93NUMPY\x01\x00");
        let hlen = u16::from_le_bytes([b[8], b[9]]) as usize;
        assert_eq!((10 + hlen) % 64, 0);
        let h = std::str::from_utf8(&b[10..10 + hlen]).unwrap();
        assert!(h.starts_with("{'descr': '<f8', 'fortran_order': False, 'shape': (2, 3), }"));
        assert!(h.ends_with('\n'));
        assert!(npy_bytes(&[4], &[0.0; 4]).windows(4).any(|w| w == b"(4,)"));
        assert_eq!(parse_npy(&b).unwrap(), (vec![2, 3], vec![0.0; 6]));
    }

    #[test]
    fn checkpoint_compatibility() {
        let c = sample(false);
        let sp = c.species.clone();
        assert!(c.check_compatible(&c.grid, &sp).is_ok());
        let mut g = c.grid.clone();
        g.dims = [2, 3, 2];
        assert!(c.check_compatible(&g, &sp).is_err());
        let mut g = c.grid.clone();
        g.spacing[1] = 2.1e-9;
        assert!(c.check_compatible(&g, &sp).is_err());
        let swapped = vec![sp[1].clone(), sp[0].clone()];
        assert!(c.check_compatible(&c.grid, &swapped).is_err());
    }

    #[test]
    fn checkpoint_read_rejects_garbage() {
        let p = tmp_path("garbage.npz");
        std::fs::write(&p, b"not a zip").unwrap();
        assert!(matches!(
            Checkpoint::read(&p),
            Err(RadkinError::Checkpoint(_))
        ));
    }
}
