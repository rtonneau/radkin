//! `[irt]` section of a configuration.
//!
//! ```toml
//! [irt]
//! t_end = 1e-6             # s (required): end of the IRT stage; the handover time in `radkin voxel`
//! t_start = 1e-12          # s (default: the earliest seed time)
//! trajectories = 16        # default 1 (`radkin irt` only)
//! random_seed = 1          # default 0 (`radkin irt` only)
//! resampling = "full"      # pair rule
//! start_contact = true     # a pair injected inside its effective radius reacts at once
//!                          # with the contact probability; false for seeds taken from a
//!                          # chemistry stage that has already started (Geant4 at 1.001 ps)
//! boundary = "none"        # "none" (free space) or "outflow" (needs a box)
//! cutoff_probability = 0.01
//! record_times = [1e-9]    # s (default: 100 log-spaced times)
//! clip = false             # drop seed molecules outside the box instead of failing
//!
//! [irt.box]                # default: the box of [voxel.grid]
//! lo_nm = [-1000, -1000, -1000]
//! hi_nm = [1000, 1000, 1000]
//!
//! [[irt.seed]]             # default: the [[voxel.seed]] entries
//! file = "molecules.csv"
//! time = 1e-12
//! events = [0]
//! offset_nm = [0, 0, 0]
//! energy_ev = 615.0
//!
//! [irt.geant4]             # as [voxel.geant4]
//! ignore = ["H2O^0"]
//!
//! [irt.species.Cl_minus]   # overrides or adds an entry of the species table
//! radius_nm = 0.18
//! charge = -1
//! diffusion = 2.0e-9       # m^2/s
//! ```

use std::path::Path;

use indexmap::IndexMap;

use crate::error::{RadkinError, Result};
use crate::irt::sample::Vec3;
use crate::irt::species::{SpeciesProps, SpeciesTable};
use crate::voxel::config::{VoxelConfig, num, seed_policy, seed_spec, triple};
use crate::voxel::seed::{SeedPolicy, SeedSpec};

const NM: f64 = 1e-9;

/// Keys of `[irt]`.
const KEYS: [&str; 14] = [
    "t_end",
    "t_start",
    "trajectories",
    "random_seed",
    "resampling",
    "start_contact",
    "boundary",
    "cutoff_probability",
    "record_times",
    "clip",
    "box",
    "geant4",
    "species",
    "seed",
];

/// Pair rule of the IRT stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Resampling {
    /// Every sampled reaction time is discarded after each reaction.
    Full,
}

/// `boundary` key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryKind {
    None,
    Outflow,
}

/// Fields of an `[irt.species.<label>]` table.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SpeciesOverride {
    /// Van der Waals radius [m].
    pub radius: Option<f64>,
    pub charge: Option<i32>,
    /// Diffusion coefficient [m^2/s].
    pub diffusion: Option<f64>,
}

/// Parsed `[irt]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct IrtConfig {
    pub t_end: f64,
    pub t_start: Option<f64>,
    /// `None` when the key is absent (default 1).
    pub trajectories: Option<usize>,
    /// `irt.random_seed`; `None` when the key is absent (default 0).
    pub seed: Option<u64>,
    pub resampling: Resampling,
    pub start_contact: bool,
    pub boundary: BoundaryKind,
    pub cutoff_probability: f64,
    pub record_times: Option<Vec<f64>>,
    /// Corners of the box [m].
    pub domain: Option<(Vec3, Vec3)>,
    pub seeds: Vec<SeedSpec>,
    pub policy: SeedPolicy,
    pub species: IndexMap<String, SpeciesOverride>,
    /// `[voxel.diffusion]` overrides [m^2/s], applied before `species`.
    pub voxel_diffusion: IndexMap<String, f64>,
}

fn boolean(t: &toml::Table, key: &str, default: bool) -> Result<bool> {
    match t.get(key) {
        None => Ok(default),
        Some(toml::Value::Boolean(b)) => Ok(*b),
        Some(_) => Err(RadkinError::Config(format!("irt.{key} must be a boolean"))),
    }
}

fn opt_num(t: &toml::Table, key: &str) -> Result<Option<f64>> {
    t.get(key)
        .map(|v| num(v, &format!("irt.{key}")))
        .transpose()
}

fn opt_count(t: &toml::Table, key: &str, min: i64) -> Result<Option<u64>> {
    match t.get(key) {
        None => Ok(None),
        Some(toml::Value::Integer(i)) if *i >= min => Ok(Some(*i as u64)),
        Some(_) => Err(RadkinError::Config(format!(
            "irt.{key} must be an integer >= {min}"
        ))),
    }
}

fn species_override(t: &toml::Table, label: &str) -> Result<SpeciesOverride> {
    let what = format!("irt.species.{label}");
    for k in t.keys() {
        if !["radius_nm", "charge", "diffusion"].contains(&k.as_str()) {
            return Err(RadkinError::Config(format!(
                "{what}: unknown key '{k}' (radius_nm, charge, diffusion)"
            )));
        }
    }
    let positive = |key: &str, allow_zero: bool| -> Result<Option<f64>> {
        match t.get(key) {
            None => Ok(None),
            Some(v) => {
                let x = num(v, &format!("{what}.{key}"))?;
                if !x.is_finite() || x < 0.0 || (x == 0.0 && !allow_zero) {
                    return Err(RadkinError::Config(format!(
                        "{what}.{key} must be finite and {} 0, got {x}",
                        if allow_zero { ">=" } else { ">" }
                    )));
                }
                Ok(Some(x))
            }
        }
    };
    let charge = match t.get("charge") {
        None => None,
        Some(toml::Value::Integer(i)) => Some(
            i32::try_from(*i)
                .map_err(|_| RadkinError::Config(format!("{what}.charge is out of range")))?,
        ),
        Some(_) => {
            return Err(RadkinError::Config(format!(
                "{what}.charge must be an integer"
            )));
        }
    };
    Ok(SpeciesOverride {
        radius: positive("radius_nm", false)?.map(|r| r * NM),
        charge,
        diffusion: positive("diffusion", true)?,
    })
}

impl IrtConfig {
    /// Read the `[irt]` section of a TOML file; `None` without one.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Option<Self>> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        Self::from_toml_str_in(&text, path.parent().unwrap_or(Path::new(".")))
            .map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
    }

    pub fn from_toml_str(text: &str) -> Result<Option<Self>> {
        Self::from_toml_str_in(text, Path::new("."))
    }

    /// Parse TOML text; relative seed files are resolved against `base`.
    pub fn from_toml_str_in(text: &str, base: &Path) -> Result<Option<Self>> {
        let root: toml::Table =
            toml::from_str(text).map_err(|e| RadkinError::Config(e.to_string()))?;
        let Some(irt) = root.get("irt") else {
            return Ok(None);
        };
        let irt = irt
            .as_table()
            .ok_or_else(|| RadkinError::Config("irt must be a table".into()))?;
        let voxel = match root.get("voxel") {
            Some(v) => Some(VoxelConfig::from_table(
                v.as_table()
                    .ok_or_else(|| RadkinError::Config("voxel must be a table".into()))?,
                base,
            )?),
            None => None,
        };
        Self::from_tables(irt, voxel.as_ref(), root.get("voxel"), base).map(Some)
    }

    fn from_tables(
        irt: &toml::Table,
        voxel: Option<&VoxelConfig>,
        voxel_raw: Option<&toml::Value>,
        base: &Path,
    ) -> Result<Self> {
        for k in irt.keys() {
            if !KEYS.contains(&k.as_str()) {
                return Err(RadkinError::Config(format!("irt: unknown key '{k}'")));
            }
        }
        let t_end = opt_num(irt, "t_end")?
            .ok_or_else(|| RadkinError::Config("irt.t_end is required".into()))?;
        let t_start = opt_num(irt, "t_start")?;
        if !t_end.is_finite() || t_start.is_some_and(|t| !(t.is_finite() && t >= 0.0 && t < t_end))
        {
            return Err(RadkinError::Config(format!(
                "irt: t_end ({t_end}) must be finite and above t_start ({t_start:?}), itself >= 0"
            )));
        }
        let resampling = match irt.get("resampling") {
            None => Resampling::Full,
            Some(toml::Value::String(s)) if s == "full" => Resampling::Full,
            Some(v) => {
                return Err(RadkinError::Config(format!(
                    "irt.resampling must be \"full\", got {v}"
                )));
            }
        };
        let boundary = match irt.get("boundary") {
            None => BoundaryKind::None,
            Some(toml::Value::String(s)) if s == "none" => BoundaryKind::None,
            Some(toml::Value::String(s)) if s == "outflow" => BoundaryKind::Outflow,
            Some(v) => {
                return Err(RadkinError::Config(format!(
                    "irt.boundary must be \"none\" or \"outflow\", got {v}"
                )));
            }
        };
        let cutoff_probability = opt_num(irt, "cutoff_probability")?.unwrap_or(0.01);
        if !(cutoff_probability > 0.0 && cutoff_probability < 1.0) {
            return Err(RadkinError::Config(format!(
                "irt.cutoff_probability must be in (0, 1), got {cutoff_probability}"
            )));
        }
        let record_times = match irt.get("record_times") {
            None => None,
            Some(v) => Some(
                v.as_array()
                    .and_then(|a| {
                        a.iter()
                            .map(|x| x.as_float().or(x.as_integer().map(|i| i as f64)))
                            .collect::<Option<Vec<f64>>>()
                    })
                    .ok_or_else(|| {
                        RadkinError::Config("irt.record_times must be a list of numbers".into())
                    })?,
            ),
        };
        let domain = match irt.get("box") {
            Some(b) => {
                let b = b
                    .as_table()
                    .ok_or_else(|| RadkinError::Config("irt.box must be a table".into()))?;
                let lo = triple(b.get("lo_nm"), None, "irt.box.lo_nm")?.map(|x| x * NM);
                let hi = triple(b.get("hi_nm"), None, "irt.box.hi_nm")?.map(|x| x * NM);
                if (0..3).any(|a| !(lo[a].is_finite() && hi[a].is_finite() && hi[a] > lo[a])) {
                    return Err(RadkinError::Config(
                        "irt.box: hi_nm must be above lo_nm on every axis".into(),
                    ));
                }
                Some((lo, hi))
            }
            None => voxel.map(|v| v.grid.bounds()),
        };
        if boundary != BoundaryKind::None && domain.is_none() {
            return Err(RadkinError::Config(
                "irt.boundary needs a box: [irt.box] or [voxel.grid]".into(),
            ));
        }
        let own_seeds = match irt.get("seed") {
            None => None,
            Some(toml::Value::Array(a)) => Some(
                a.iter()
                    .enumerate()
                    .map(|(i, v)| {
                        let t = v.as_table().ok_or_else(|| {
                            RadkinError::Config("irt.seed entries must be tables".into())
                        })?;
                        seed_spec(t, base, &format!("irt.seed[{i}]"))
                    })
                    .collect::<Result<Vec<_>>>()?,
            ),
            Some(_) => {
                return Err(RadkinError::Config(
                    "irt.seed must be an array of tables ([[irt.seed]])".into(),
                ));
            }
        };
        let seeds = match (own_seeds, voxel) {
            (Some(s), _) => s,
            (None, Some(v)) => v.seeds.clone(),
            (None, None) => Vec::new(),
        };
        // [irt.geant4] and irt.clip, else those of [voxel].
        let policy = if irt.contains_key("geant4") || irt.contains_key("clip") || voxel.is_none() {
            seed_policy(irt, "irt")?
        } else {
            voxel.map(|v| v.policy.clone()).unwrap_or_default()
        };
        let mut species = IndexMap::new();
        if let Some(sp) = irt.get("species") {
            let sp = sp
                .as_table()
                .ok_or_else(|| RadkinError::Config("irt.species must be a table".into()))?;
            for (label, v) in sp {
                let t = v.as_table().ok_or_else(|| {
                    RadkinError::Config(format!("irt.species.{label} must be a table"))
                })?;
                species.insert(label.clone(), species_override(t, label)?);
            }
        }
        // Only the coefficients the config sets, not the defaults of [voxel].
        let mut voxel_diffusion = IndexMap::new();
        if let (Some(v), Some(raw)) = (voxel, voxel_raw)
            && let Some(d) = raw.get("diffusion").and_then(|d| d.as_table())
        {
            for k in d.keys() {
                voxel_diffusion.insert(k.clone(), v.diffusion_of(k));
            }
        }
        Ok(Self {
            t_end,
            t_start,
            trajectories: opt_count(irt, "trajectories", 1)?.map(|n| n as usize),
            seed: opt_count(irt, "random_seed", 0)?,
            resampling,
            start_contact: boolean(irt, "start_contact", true)?,
            boundary,
            cutoff_probability,
            record_times,
            domain,
            seeds,
            policy,
            species,
            voxel_diffusion,
        })
    }

    /// The species table of the run: the built-in water table, the
    /// `[voxel.diffusion]` overrides, then `[irt.species]`.
    ///
    /// A species the built-in table does not hold needs a radius and a
    /// charge to get an entry; with a diffusion coefficient only, it is left
    /// out of the table and its coefficient is returned in the second map.
    pub fn species_table(&self) -> Result<(SpeciesTable, IndexMap<String, f64>)> {
        let mut table = SpeciesTable::water();
        let mut extra: IndexMap<String, f64> = IndexMap::new();
        for (label, d) in &self.voxel_diffusion {
            match table.get(label).copied() {
                Some(p) => table.set(label, SpeciesProps { diffusion: *d, ..p }),
                None => {
                    extra.insert(label.clone(), *d);
                }
            }
        }
        for (label, o) in &self.species {
            match (table.get(label).copied(), o.radius, o.charge) {
                (Some(p), _, _) => table.set(
                    label,
                    SpeciesProps {
                        radius: o.radius.unwrap_or(p.radius),
                        charge: o.charge.unwrap_or(p.charge),
                        diffusion: o.diffusion.unwrap_or(p.diffusion),
                    },
                ),
                (None, Some(radius), Some(charge)) => {
                    let diffusion = o
                        .diffusion
                        .or_else(|| extra.get(label).copied())
                        .ok_or_else(|| {
                            RadkinError::Config(format!(
                                "irt.species.{label}: diffusion is required for a species the built-in table does not hold"
                            ))
                        })?;
                    extra.shift_remove(label);
                    table.set(
                        label,
                        SpeciesProps {
                            radius,
                            charge,
                            diffusion,
                        },
                    );
                }
                (None, None, None) => {
                    if let Some(d) = o.diffusion {
                        extra.insert(label.clone(), d);
                    }
                }
                (None, _, _) => {
                    return Err(RadkinError::Config(format!(
                        "irt.species.{label}: radius_nm and charge must be given together for a species the built-in table does not hold"
                    )));
                }
            }
        }
        Ok((table, extra))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn absent_section() {
        assert_eq!(
            IrtConfig::from_toml_str("[parameters]\npH = 7").unwrap(),
            None
        );
    }

    #[test]
    fn defaults() {
        let c = IrtConfig::from_toml_str("[irt]\nt_end = 1e-6")
            .unwrap()
            .unwrap();
        assert_eq!(c.t_end, 1e-6);
        assert_eq!(c.t_start, None);
        assert_eq!((c.trajectories, c.seed), (None, None));
        assert_eq!(c.resampling, Resampling::Full);
        assert!(c.start_contact);
        assert_eq!(c.boundary, BoundaryKind::None);
        assert_eq!(c.cutoff_probability, 0.01);
        assert_eq!(c.record_times, None);
        assert_eq!(c.domain, None);
        assert!(c.seeds.is_empty());
        assert_eq!(c.policy.aliases["e_aq^-1"], "e_aq");
    }

    #[test]
    fn parses_every_key() {
        let c = IrtConfig::from_toml_str_in(
            r#"
            [irt]
            t_end = 1e-5
            t_start = 1.001e-12
            trajectories = 8
            random_seed = 3
            resampling = "full"
            start_contact = false
            boundary = "outflow"
            cutoff_probability = 0.001
            record_times = [1e-9, 1e-6]
            clip = true
            [irt.box]
            lo_nm = -1000
            hi_nm = [1000, 1000, 2000]
            [irt.geant4]
            ignore = ["H2O^0"]
            [[irt.seed]]
            file = "m.csv"
            time = 1.001e-12
            events = [4]
            energy_ev = 500
            [irt.species.e_aq]
            radius_nm = 0.4
            [irt.species.Cl_minus]
            radius_nm = 0.18
            charge = -1
            diffusion = 2e-9
            "#,
            Path::new("base"),
        )
        .unwrap()
        .unwrap();
        assert_eq!(c.t_start, Some(1.001e-12));
        assert_eq!((c.trajectories, c.seed), (Some(8), Some(3)));
        assert!(!c.start_contact);
        assert_eq!(c.boundary, BoundaryKind::Outflow);
        assert_eq!(c.record_times, Some(vec![1e-9, 1e-6]));
        let (lo, hi) = c.domain.unwrap();
        assert!((lo[0] + 1e-6).abs() < 1e-18 && (hi[2] - 2e-6).abs() < 1e-18);
        assert!(c.policy.clip);
        assert_eq!(c.policy.ignore, vec!["H2O^0"]);
        assert_eq!(c.seeds.len(), 1);
        assert_eq!(c.seeds[0].file, Path::new("base").join("m.csv"));
        assert_eq!(c.seeds[0].events, Some(vec![4]));
        let (table, extra) = c.species_table().unwrap();
        assert!((table.get("e_aq").unwrap().radius - 0.4e-9).abs() < 1e-20);
        assert_eq!(table.get("e_aq").unwrap().diffusion, 4.9e-9);
        assert_eq!(table.get("Cl_minus").unwrap().charge, -1);
        assert!(extra.is_empty());
    }

    #[test]
    fn falls_back_to_the_voxel_section() {
        let c = IrtConfig::from_toml_str(
            r#"
            [irt]
            t_end = 1e-6
            boundary = "outflow"
            [voxel]
            clip = true
            [voxel.grid]
            dims = [10, 10, 20]
            spacing_nm = 100
            origin_nm = [-500, -500, -1000]
            [voxel.diffusion]
            e_aq = 5e-9
            HO3 = 1e-9
            [[voxel.seed]]
            file = "m.csv"
            time = 1e-12
            "#,
        )
        .unwrap()
        .unwrap();
        let (lo, hi) = c.domain.unwrap();
        assert!((lo[2] + 1e-6).abs() < 1e-18 && (hi[0] - 0.5e-6).abs() < 1e-18);
        assert_eq!(c.seeds.len(), 1);
        assert!(c.policy.clip);
        let (table, extra) = c.species_table().unwrap();
        assert_eq!(table.get("e_aq").unwrap().diffusion, 5e-9);
        assert_eq!(table.get("OH_r").unwrap().diffusion, 2.2e-9);
        // not in the built-in table: a diffusion coefficient only
        assert!(table.get("HO3").is_none());
        assert_eq!(extra["HO3"], 1e-9);
    }

    #[test]
    fn rejects_bad_values() {
        for (bad, msg) in [
            ("[irt]", "t_end is required"),
            ("[irt]\nt_end = 1e-6\nt_start = 1e-5", "above t_start"),
            ("[irt]\nt_end = 1e-6\nresampling = \"lazy\"", "resampling"),
            ("[irt]\nt_end = 1e-6\nboundary = \"outflow\"", "needs a box"),
            ("[irt]\nt_end = 1e-6\nboundary = \"periodic\"", "boundary"),
            ("[irt]\nt_end = 1e-6\ntrajectories = 0", "trajectories"),
            (
                "[irt]\nt_end = 1e-6\ncutoff_probability = 2",
                "cutoff_probability",
            ),
            ("[irt]\nt_end = 1e-6\nhandover_time = 1", "unknown key"),
            (
                "[irt]\nt_end = 1e-6\n[irt.box]\nlo_nm = 1\nhi_nm = 0",
                "hi_nm must be above",
            ),
            (
                "[irt]\nt_end = 1e-6\n[irt.species.e_aq]\nradius = 1",
                "unknown key",
            ),
        ] {
            let e = IrtConfig::from_toml_str(bad).unwrap_err().to_string();
            assert!(e.contains(msg), "{bad}: {e}");
        }
        let c = IrtConfig::from_toml_str("[irt]\nt_end = 1\n[irt.species.X_r]\nradius_nm = 0.2")
            .unwrap()
            .unwrap();
        let e = c.species_table().unwrap_err().to_string();
        assert!(e.contains("given together"), "{e}");
    }
}
