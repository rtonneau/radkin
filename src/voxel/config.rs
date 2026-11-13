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
use crate::voxel::seed::{SeedPolicy, SeedSpec};

/// Default diffusion coefficients [m^2/s] of the radkin water species,
/// from Geant4-DNA 11.4.3 option3: the configuration-level values set in
/// `G4ChemDissociationChannels_option1::ConstructMolecule` where it sets one
/// (`H3Op`, `°OH`, `OHm`, `H2`, `H2O2`, `HO2m`, `O2m`, `Om`, `O3m`), otherwise
/// the molecule definition's (`G4Electron_aq.cc`, `G4Hydrogen.cc`, `G4HO2.cc`,
/// `G4O2.cc`, `G4O3.cc`, `G4Oxygen.cc`). Geant4-DNA has no HO3 (a species of the
/// Tonneau 2025 chemistry): it takes O3's 2.0e-9, an assumption. Species absent here
/// and from `[voxel.diffusion]` do not diffuse.
pub const DEFAULT_DIFFUSION: [(&str, f64); 16] = [
    ("e_aq", 4.9e-9),
    ("OH_r", 2.2e-9),
    ("H3O_plus", 9.46e-9),
    ("H_r", 7.0e-9),
    ("H2", 4.8e-9),
    ("H2O2", 2.3e-9),
    ("OH_minus", 5.3e-9),
    ("HO2", 2.3e-9),
    ("HO2_minus", 1.4e-9),
    ("O2", 2.4e-9),
    ("O2_minus", 1.75e-9),
    ("O3", 2.0e-9),
    ("O3_minus", 2.0e-9),
    ("O_minus", 2.0e-9),
    ("O_atom", 2.0e-9),
    ("HO3", 2.0e-9), // not in Geant4-DNA: O3's value
];

const NM: f64 = 1e-9;

/// Parsed `[voxel]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelConfig {
    pub grid: VoxelGrid,
    /// Diffusion coefficient [m^2/s] per species label (defaults + overrides).
    pub diffusion: IndexMap<String, f64>,
    /// `[[voxel.seed]]` entries (files resolved against the config folder).
    pub seeds: Vec<SeedSpec>,
    /// `[voxel.geant4]` aliases and ignore list, plus `voxel.clip`.
    pub policy: SeedPolicy,
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

fn strings(v: &toml::Value, what: &str) -> Result<Vec<String>> {
    v.as_array()
        .and_then(|a| a.iter().map(|x| x.as_str().map(String::from)).collect())
        .ok_or_else(|| RadkinError::Config(format!("{what} must be a list of strings")))
}

fn seed_spec(t: &toml::Table, base: &Path, i: usize) -> Result<SeedSpec> {
    let what = format!("voxel.seed[{i}]");
    let file = t
        .get("file")
        .and_then(|v| v.as_str())
        .ok_or_else(|| RadkinError::Config(format!("{what}.file (string) is required")))?;
    let file = Path::new(file);
    let file = if file.is_absolute() {
        file.to_path_buf()
    } else {
        base.join(file)
    };
    let time = num(
        t.get("time")
            .ok_or_else(|| RadkinError::Config(format!("{what}.time is required")))?,
        &format!("{what}.time"),
    )?;
    if !(time.is_finite() && time >= 0.0) {
        return Err(RadkinError::Config(format!("{what}.time must be >= 0")));
    }
    let events = t
        .get("events")
        .map(|v| {
            v.as_array()
                .and_then(|a| a.iter().map(|x| x.as_integer()).collect::<Option<Vec<_>>>())
                .ok_or_else(|| RadkinError::Config(format!("{what}.events must be integers")))
        })
        .transpose()?;
    let offset =
        triple(t.get("offset_nm"), Some(0.0), &format!("{what}.offset_nm"))?.map(|x| x * NM);
    let energy_ev = t
        .get("energy_ev")
        .map(|v| num(v, &format!("{what}.energy_ev")))
        .transpose()?;
    Ok(SeedSpec {
        file,
        time,
        events,
        offset,
        energy_ev,
    })
}

impl VoxelConfig {
    /// Read the `[voxel]` section of a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        Self::from_toml_str_in(&text, path.parent().unwrap_or(Path::new(".")))
            .map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
    }

    pub fn from_toml_str(text: &str) -> Result<Self> {
        Self::from_toml_str_in(text, Path::new("."))
    }

    /// Parse TOML text; relative seed files are resolved against `base`.
    pub fn from_toml_str_in(text: &str, base: &Path) -> Result<Self> {
        let root: toml::Table =
            toml::from_str(text).map_err(|e| RadkinError::Config(e.to_string()))?;
        let voxel = root
            .get("voxel")
            .and_then(|v| v.as_table())
            .ok_or_else(|| RadkinError::Config("missing [voxel] section".into()))?;
        Self::from_table(voxel, base)
    }

    pub fn from_table(voxel: &toml::Table, base: &Path) -> Result<Self> {
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
        let seeds = match voxel.get("seed") {
            None => Vec::new(),
            Some(toml::Value::Array(a)) => a
                .iter()
                .enumerate()
                .map(|(i, v)| {
                    let t = v.as_table().ok_or_else(|| {
                        RadkinError::Config("voxel.seed entries must be tables".into())
                    })?;
                    seed_spec(t, base, i)
                })
                .collect::<Result<_>>()?,
            Some(_) => {
                return Err(RadkinError::Config(
                    "voxel.seed must be an array of tables ([[voxel.seed]])".into(),
                ));
            }
        };
        let mut policy = SeedPolicy {
            clip: voxel
                .get("clip")
                .map(|v| {
                    v.as_bool()
                        .ok_or_else(|| RadkinError::Config("voxel.clip must be a boolean".into()))
                })
                .transpose()?
                .unwrap_or(false),
            ..SeedPolicy::default()
        };
        if let Some(g4) = voxel.get("geant4") {
            let g4 = g4
                .as_table()
                .ok_or_else(|| RadkinError::Config("voxel.geant4 must be a table".into()))?;
            if let Some(a) = g4.get("aliases") {
                let a = a.as_table().ok_or_else(|| {
                    RadkinError::Config("voxel.geant4.aliases must be a table".into())
                })?;
                for (k, v) in a {
                    let v = v.as_str().ok_or_else(|| {
                        RadkinError::Config(format!("voxel.geant4.aliases.{k} must be a string"))
                    })?;
                    policy.aliases.insert(k.clone(), v.to_string());
                }
            }
            if let Some(i) = g4.get("ignore") {
                policy.ignore = strings(i, "voxel.geant4.ignore")?;
            }
        }
        Ok(Self {
            grid,
            diffusion,
            seeds,
            policy,
        })
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
        assert_eq!(c.diffusion_of("OH_r"), 2.2e-9);
        assert_eq!(c.diffusion_of("HO3"), 2.0e-9);
        assert_eq!(c.diffusion_of("H2O"), 0.0);
    }

    #[test]
    fn parses_seeds_and_policy() {
        let c = VoxelConfig::from_toml_str_in(
            r#"
            [voxel]
            clip = true
            [voxel.grid]
            dims = [1, 1, 1]
            spacing_nm = 1
            [voxel.geant4]
            aliases = { "O2^0" = "O2" }
            ignore = ["H2O^0"]
            [[voxel.seed]]
            file = "track.csv"
            time = 1e-12
            events = [0, 2]
            offset_nm = [1, 0, 0]
            energy_ev = 1000.0
            [[voxel.seed]]
            file = "track.csv"
            time = 1e-6
            "#,
            Path::new("base"),
        )
        .unwrap();
        assert!(c.policy.clip);
        assert_eq!(c.policy.aliases["O2^0"], "O2");
        assert_eq!(c.policy.aliases["e_aq^-1"], "e_aq");
        assert_eq!(c.policy.ignore, vec!["H2O^0"]);
        assert_eq!(c.seeds.len(), 2);
        assert_eq!(c.seeds[0].file, Path::new("base").join("track.csv"));
        assert_eq!(c.seeds[0].events, Some(vec![0, 2]));
        assert_eq!(c.seeds[0].energy_ev, Some(1000.0));
        assert_eq!(c.seeds[1].time, 1e-6);
        assert_eq!(c.seeds[1].events, None);
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
