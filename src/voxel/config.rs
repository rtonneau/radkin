//! `[voxel]` section of a run TOML.
//!
//! ```toml
//! [voxel.grid]
//! dims = [20, 20, 20]
//! spacing_nm = [5.0, 5.0, 5.0]   # or one number
//! origin_nm = [-50.0, -50.0, -50.0]  # default 0
//!
//! [voxel.diffusion]   # m^2/s, overrides DEFAULT_DIFFUSION
//! O2 = 2.4e-9
//! ```

use std::path::Path;

use indexmap::IndexMap;

use crate::error::{RadkinError, Result};
use crate::voxel::grid::VoxelGrid;

/// Default diffusion coefficients [m^2/s] of the radkin water species,
/// copied from Geant4-DNA 11.4.3 (`G4Electron_aq.cc`, `G4OH.cc`, `G4H3O.cc`,
/// `G4Hydrogen.cc`, `G4H2.cc`, `G4H2O2.cc`; `OHm` from
/// `G4ChemDissociationChannels.cc`). Species absent here and from
/// `[voxel.diffusion]` do not diffuse.
pub const DEFAULT_DIFFUSION: [(&str, f64); 7] = [
    ("e_aq", 4.9e-9),
    ("OH_r", 2.8e-9),
    ("H3O_plus", 9.0e-9),
    ("H_r", 7.0e-9),
    ("H2", 5.0e-9),
    ("H2O2", 1.4e-9),
    ("OH_minus", 5.0e-9),
];

const NM: f64 = 1e-9;

/// Parsed `[voxel]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelConfig {
    pub grid: VoxelGrid,
    /// Diffusion coefficient [m^2/s] per species label (defaults + overrides).
    pub diffusion: IndexMap<String, f64>,
}

fn num(v: &toml::Value, what: &str) -> Result<f64> {
    match v {
        toml::Value::Float(f) => Ok(*f),
        toml::Value::Integer(i) => Ok(*i as f64),
        _ => Err(RadkinError::Config(format!("{what} must be a number"))),
    }
}

fn triple(v: Option<&toml::Value>, default: Option<f64>, what: &str) -> Result<[f64; 3]> {
    match v {
        None => default
            .map(|d| [d; 3])
            .ok_or_else(|| RadkinError::Config(format!("{what} is required"))),
        Some(toml::Value::Array(a)) if a.len() == 3 => {
            Ok([num(&a[0], what)?, num(&a[1], what)?, num(&a[2], what)?])
        }
        Some(x @ (toml::Value::Float(_) | toml::Value::Integer(_))) => Ok([num(x, what)?; 3]),
        Some(_) => Err(RadkinError::Config(format!(
            "{what} must be a number or a list of 3 numbers"
        ))),
    }
}

impl VoxelConfig {
    /// Read the `[voxel]` section of a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        Self::from_toml_str(&text)
            .map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
    }

    pub fn from_toml_str(text: &str) -> Result<Self> {
        let root: toml::Table =
            toml::from_str(text).map_err(|e| RadkinError::Config(e.to_string()))?;
        let voxel = root
            .get("voxel")
            .and_then(|v| v.as_table())
            .ok_or_else(|| RadkinError::Config("missing [voxel] section".into()))?;
        Self::from_table(voxel)
    }

    pub fn from_table(voxel: &toml::Table) -> Result<Self> {
        let g = voxel
            .get("grid")
            .and_then(|v| v.as_table())
            .ok_or_else(|| RadkinError::Config("missing [voxel.grid] section".into()))?;
        let dims = match g.get("dims") {
            Some(toml::Value::Array(a)) if a.len() == 3 => {
                let mut d = [0usize; 3];
                for (o, v) in d.iter_mut().zip(a) {
                    *o = v.as_integer().filter(|i| *i >= 0).ok_or_else(|| {
                        RadkinError::Config("voxel.grid.dims must be non-negative integers".into())
                    })? as usize;
                }
                d
            }
            _ => {
                return Err(RadkinError::Config(
                    "voxel.grid.dims must be a list of 3 integers".into(),
                ));
            }
        };
        let spacing = triple(g.get("spacing_nm"), None, "voxel.grid.spacing_nm")?.map(|x| x * NM);
        let origin = triple(g.get("origin_nm"), Some(0.0), "voxel.grid.origin_nm")?.map(|x| x * NM);
        let grid = VoxelGrid::new(dims, spacing, origin)?;

        let mut diffusion: IndexMap<String, f64> = DEFAULT_DIFFUSION
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        if let Some(t) = voxel.get("diffusion") {
            let t = t
                .as_table()
                .ok_or_else(|| RadkinError::Config("voxel.diffusion must be a table".into()))?;
            for (k, v) in t {
                let d = num(v, &format!("voxel.diffusion.{k}"))?;
                if !(d.is_finite() && d >= 0.0) {
                    return Err(RadkinError::Config(format!(
                        "voxel.diffusion.{k} must be finite and >= 0, got {d}"
                    )));
                }
                diffusion.insert(k.clone(), d);
            }
        }
        Ok(Self { grid, diffusion })
    }

    /// Diffusion coefficient [m^2/s] of `label` (0 = immobile).
    pub fn diffusion_of(&self, label: &str) -> f64 {
        self.diffusion.get(label).copied().unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_grid_and_diffusion() {
        let c = VoxelConfig::from_toml_str(
            r#"
            [voxel.grid]
            dims = [4, 5, 6]
            spacing_nm = 2
            origin_nm = [-4, -5, -6]
            [voxel.diffusion]
            O2 = 2.4e-9
            e_aq = 5e-9
            "#,
        )
        .unwrap();
        assert_eq!(c.grid.dims, [4, 5, 6]);
        let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-24);
        assert!(close(c.grid.spacing, [2e-9; 3]));
        assert!(close(c.grid.origin, [-4e-9, -5e-9, -6e-9]));
        assert_eq!(c.diffusion_of("O2"), 2.4e-9);
        assert_eq!(c.diffusion_of("e_aq"), 5e-9);
        assert_eq!(c.diffusion_of("OH_r"), 2.8e-9);
        assert_eq!(c.diffusion_of("H2O"), 0.0);
    }

    #[test]
    fn rejects_invalid() {
        for bad in [
            "[voxel.grid]\ndims = [0, 1, 1]\nspacing_nm = 1",
            "[voxel.grid]\ndims = [1, 1, 1]\nspacing_nm = 0",
            "[voxel.grid]\ndims = [1, 1]\nspacing_nm = 1",
            "[voxel.grid]\ndims = [1, 1, 1]",
            "[voxel.grid]\ndims = [1, 1, 1]\nspacing_nm = 1\n[voxel.diffusion]\nO2 = -1.0",
            "[parameters]\npH = 7",
        ] {
            assert!(VoxelConfig::from_toml_str(bad).is_err(), "{bad}");
        }
    }
}
