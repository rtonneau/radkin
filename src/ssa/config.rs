//! Solver selection and the `[ssa]` table of a run TOML.
//!
//! ```toml
//! solver = "ssa"        # "lsoda" (default) or "ssa"; before any [table]
//!
//! [ssa]
//! volume = 1e-15        # L, unit-cell volume (voxel runs use the voxel volume)
//! trajectories = 100    # ensemble size
//! seed = 0              # trajectory i draws from stream i of this seed
//! tau_leap = false      # opt-in tau-leaping
//! tau_epsilon = 0.03    # accuracy of a leap (smaller: closer to the exact SSA)
//! rate_correction = "hellander"  # voxel runs: "hellander" (default) or "none"
//! background_run = true # voxel runs: subtract a seedless run (default true)
//!
//! [ssa.species_levels]  # voxel runs: grid level of a species (default 0, fine)
//! O2 = 2                # O2 hops and reacts on 4x4x4 blocks of fine voxels
//! ```
//!
//! `rate_correction` scales the bimolecular propensities of a voxel run for
//! the voxel size (mesh correction, `docs/ssa.md`). Unit-cell runs have no
//! grid: an explicit `"hellander"` there is an error.
//!
//! `background_run` runs a second, seedless ensemble on the same grid and
//! settings, with a seed derived from `seed`, and subtracts its counts so
//! that the G-values of `yields.csv` are background-net (`N_` stays gross). A
//! unit-cell run has no seeds: `true` there is an error, `false` is accepted.
//!
//! `species_levels` puts a species on a coarser level of the 2^L voxel
//! hierarchy (multi-grid RDME, ADR 0005). It needs a stochastic voxel run.
//!
//! Both are read from the run file itself, not from its includes.

use std::path::Path;

use crate::error::{RadkinError, Result};

/// Which solver integrates the chemistry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SolverKind {
    /// Deterministic rate equations with LSODA (rsoda).
    #[default]
    Lsoda,
    /// Stochastic simulation on molecule counts.
    Ssa,
}

impl SolverKind {
    /// `"lsoda"` or `"ssa"`.
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "lsoda" => Ok(Self::Lsoda),
            "ssa" => Ok(Self::Ssa),
            other => Err(RadkinError::Config(format!(
                "solver must be \"lsoda\" or \"ssa\", got \"{other}\""
            ))),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::Lsoda => "lsoda",
            Self::Ssa => "ssa",
        }
    }
}

/// Default ensemble size.
pub const DEFAULT_TRAJECTORIES: usize = 100;

/// Default bound on the relative change of a propensity over a leap
/// (Cao, Gillespie and Petzold 2006).
pub const DEFAULT_TAU_EPSILON: f64 = 0.03;

/// Correction of the bimolecular propensities of the RDME for the voxel
/// size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RateCorrection {
    /// Macroscopic rate constants in every voxel.
    None,
    /// `c / (1 - h*/h)` on a level of spacing `h` above the critical spacing
    /// `h*` (Hellander and Hellander 2020, eq. 6).
    Hellander,
}

impl RateCorrection {
    /// `"none"` or `"hellander"`.
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "none" => Ok(Self::None),
            "hellander" => Ok(Self::Hellander),
            other => Err(RadkinError::Config(format!(
                "ssa.rate_correction must be \"hellander\" or \"none\", got \"{other}\""
            ))),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Hellander => "hellander",
        }
    }
}

/// `[ssa]` table.
#[derive(Debug, Clone, PartialEq)]
pub struct SsaConfig {
    /// Unit-cell volume [L]; required by unit-cell runs.
    pub volume: Option<f64>,
    /// Number of trajectories of the ensemble (>= 1).
    pub trajectories: usize,
    /// Seed of the generator; trajectory `i` uses its stream `i`.
    pub seed: u64,
    /// Tau-leaping instead of the exact SSA.
    pub tau_leap: bool,
    /// Bound on the relative change of a propensity over a leap, in (0, 1).
    /// The error of tau-leaping on the mean is proportional to it.
    pub tau_epsilon: f64,
    /// Mesh correction of a voxel run; `None`: not set, the default
    /// (Hellander where the grid allows it).
    pub rate_correction: Option<RateCorrection>,
    /// Seedless background run of a voxel run; `None`: not set, the default
    /// (on).
    pub background_run: Option<bool>,
    /// `[ssa.species_levels]`: `(species label, level)` in file order.
    pub species_levels: Vec<(String, u32)>,
}

/// Largest accepted species level: 2^16 fine voxels per coarse voxel along an axis.
pub const MAX_SPECIES_LEVEL: u32 = 16;

impl Default for SsaConfig {
    fn default() -> Self {
        Self {
            volume: None,
            trajectories: DEFAULT_TRAJECTORIES,
            seed: 0,
            tau_leap: false,
            tau_epsilon: DEFAULT_TAU_EPSILON,
            rate_correction: None,
            background_run: None,
            species_levels: Vec::new(),
        }
    }
}

impl SsaConfig {
    pub fn from_table(t: &toml::Table) -> Result<Self> {
        let mut out = Self::default();
        for (key, value) in t {
            match key.as_str() {
                "volume" => {
                    let v = match value {
                        toml::Value::Float(f) => *f,
                        toml::Value::Integer(i) => *i as f64,
                        _ => f64::NAN,
                    };
                    if !(v.is_finite() && v > 0.0) {
                        return Err(RadkinError::Config(format!(
                            "ssa.volume must be a number > 0 (litres), got {value}"
                        )));
                    }
                    out.volume = Some(v);
                }
                "trajectories" => match value.as_integer() {
                    Some(n) if n >= 1 => out.trajectories = n as usize,
                    _ => {
                        return Err(RadkinError::Config(format!(
                            "ssa.trajectories must be an integer >= 1, got {value}"
                        )));
                    }
                },
                "seed" => match value.as_integer() {
                    Some(n) if n >= 0 => out.seed = n as u64,
                    _ => {
                        return Err(RadkinError::Config(format!(
                            "ssa.seed must be an integer >= 0, got {value}"
                        )));
                    }
                },
                "tau_leap" => match value.as_bool() {
                    Some(b) => out.tau_leap = b,
                    None => {
                        return Err(RadkinError::Config(format!(
                            "ssa.tau_leap must be a boolean, got {value}"
                        )));
                    }
                },
                "tau_epsilon" => match value.as_float() {
                    Some(e) if e > 0.0 && e < 1.0 => out.tau_epsilon = e,
                    _ => {
                        return Err(RadkinError::Config(format!(
                            "ssa.tau_epsilon must be a number in (0, 1), got {value}"
                        )));
                    }
                },
                "rate_correction" => match value.as_str() {
                    Some(name) => out.rate_correction = Some(RateCorrection::parse(name)?),
                    None => {
                        return Err(RadkinError::Config(format!(
                            "ssa.rate_correction must be \"hellander\" or \"none\", got {value}"
                        )));
                    }
                },
                "background_run" => match value.as_bool() {
                    Some(b) => out.background_run = Some(b),
                    None => {
                        return Err(RadkinError::Config(format!(
                            "ssa.background_run must be a boolean, got {value}"
                        )));
                    }
                },
                "species_levels" => out.species_levels = parse_species_levels(value)?,
                other => {
                    return Err(RadkinError::Config(format!(
                        "unknown key ssa.{other} (known: volume, trajectories, seed, tau_leap, \
                         tau_epsilon, rate_correction, background_run, species_levels)"
                    )));
                }
            }
        }
        Ok(out)
    }

    /// The unit-cell volume [L], an error if `volume` is absent or if the
    /// run asks for a mesh correction or species levels (a unit cell has no
    /// grid).
    pub fn require_volume(&self) -> Result<f64> {
        self.refuse_species_levels("a unit-cell run")?;
        if self.background_run == Some(true) {
            return Err(RadkinError::Config(
                "ssa.background_run needs a voxel run: a unit cell has no seeds".into(),
            ));
        }
        if self.rate_correction == Some(RateCorrection::Hellander) {
            return Err(RadkinError::Config(
                "ssa.rate_correction = \"hellander\" needs a voxel run: a unit cell has no grid"
                    .into(),
            ));
        }
        self.volume.ok_or_else(|| {
            RadkinError::Config(
                "ssa.volume (litres) is required by a unit-cell run with solver = \"ssa\"".into(),
            )
        })
    }
}

impl SsaConfig {
    /// An error naming `run` if `[ssa.species_levels]` is set: only the
    /// stochastic voxel run keeps species on their own grid level.
    pub fn refuse_species_levels(&self, run: &str) -> Result<()> {
        match self.species_levels.first() {
            None => Ok(()),
            Some((name, _)) => Err(RadkinError::Config(format!(
                "ssa.species_levels ({name}) needs a voxel run with solver = \"ssa\", not {run}"
            ))),
        }
    }

    /// The level of every dynamic species, in state order (0 where not
    /// set). An unknown label is an error.
    pub fn resolve_species_levels(&self, species: &[String]) -> Result<Vec<u32>> {
        let mut out = vec![0; species.len()];
        for (name, level) in &self.species_levels {
            let Some(s) = species.iter().position(|x| x == name) else {
                return Err(RadkinError::Config(format!(
                    "ssa.species_levels.{name}: no dynamic species of that name (known: {})",
                    species.join(", ")
                )));
            };
            out[s] = *level;
        }
        Ok(out)
    }
}

/// `[ssa.species_levels]`: a table of integer levels in `0..=MAX_SPECIES_LEVEL`.
fn parse_species_levels(value: &toml::Value) -> Result<Vec<(String, u32)>> {
    let Some(t) = value.as_table() else {
        return Err(RadkinError::Config(format!(
            "ssa.species_levels must be a table of species = level, got {value}"
        )));
    };
    t.iter()
        .map(|(name, v)| match v.as_integer() {
            Some(l) if (0..=MAX_SPECIES_LEVEL as i64).contains(&l) => Ok((name.clone(), l as u32)),
            _ => Err(RadkinError::Config(format!(
                "ssa.species_levels.{name} must be an integer in 0..={MAX_SPECIES_LEVEL}, got {v}"
            ))),
        })
        .collect()
}

/// Levels of the voxel hierarchy a stochastic run needs: those of the
/// global coarsening and those of the species.
pub fn hierarchy_depth(max_coarsen_levels: u32, species_levels: &[u32]) -> u32 {
    species_levels
        .iter()
        .copied()
        .fold(max_coarsen_levels, u32::max)
}

/// Suggested level of a background species of diffusion coefficient `d`,
/// `d_slowest` being the smallest nonzero coefficient of the seeded
/// species: `floor(log2 sqrt(d / d_slowest))`, between 0 and `depth`. The
/// spacing ratio `sqrt(D_A / D_B)` equalizes the hop rates per molecule
/// (Erban and Winkelmann 2024). `None` unless both are positive.
pub fn suggested_level(d: f64, d_slowest: f64, depth: u32) -> Option<u32> {
    if !(d > 0.0 && d_slowest > 0.0 && d.is_finite() && d_slowest.is_finite()) {
        return None;
    }
    let l = (0.5 * (d / d_slowest).log2()).floor();
    Some(l.clamp(0.0, depth as f64) as u32)
}

/// `solver` and `[ssa]` of a run file.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SolverConfig {
    pub solver: SolverKind,
    pub ssa: SsaConfig,
}

impl SolverConfig {
    pub fn from_table(root: &toml::Table) -> Result<Self> {
        let solver = match root.get("solver") {
            None => SolverKind::Lsoda,
            Some(toml::Value::String(s)) => SolverKind::parse(s)?,
            Some(v) => {
                return Err(RadkinError::Config(format!(
                    "solver must be \"lsoda\" or \"ssa\", got {v}"
                )));
            }
        };
        let ssa = match root.get("ssa") {
            None => SsaConfig::default(),
            Some(toml::Value::Table(t)) => SsaConfig::from_table(t)?,
            Some(_) => return Err(RadkinError::Config("[ssa] must be a table".into())),
        };
        Ok(Self { solver, ssa })
    }

    pub fn from_toml_str(text: &str) -> Result<Self> {
        let root: toml::Table =
            toml::from_str(text).map_err(|e| RadkinError::Config(e.to_string()))?;
        Self::from_table(&root)
    }

    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        Self::from_toml_str(&text)
            .map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_to_lsoda() {
        let c = SolverConfig::from_toml_str("[parameters]\npH = 7").unwrap();
        assert_eq!(c.solver, SolverKind::Lsoda);
        assert_eq!(c.ssa, SsaConfig::default());
        assert!(c.ssa.require_volume().is_err());
    }

    #[test]
    fn parses_solver_and_ssa_table() {
        let c = SolverConfig::from_toml_str(
            "solver = \"ssa\"\n[ssa]\nvolume = 2e-15\ntrajectories = 8\nseed = 42\ntau_leap = true
tau_epsilon = 0.01\nrate_correction = \"none\"",
        )
        .unwrap();
        assert_eq!(c.solver, SolverKind::Ssa);
        assert_eq!(
            c.ssa,
            SsaConfig {
                volume: Some(2e-15),
                trajectories: 8,
                seed: 42,
                tau_leap: true,
                tau_epsilon: 0.01,
                rate_correction: Some(RateCorrection::None),
                background_run: None,
                species_levels: Vec::new(),
            }
        );
        assert_eq!(c.ssa.require_volume().unwrap(), 2e-15);
    }

    #[test]
    fn rejects_invalid_values() {
        for bad in [
            "solver = \"gillespie\"",
            "solver = 3",
            "ssa = 1",
            "[ssa]\nvolume = 0.0",
            "[ssa]\nvolume = \"big\"",
            "[ssa]\ntrajectories = 0",
            "[ssa]\nseed = -1",
            "[ssa]\ntau_leap = 1",
            "[ssa]\nvolumes = 1e-15",
            "[ssa]\nrate_correction = \"fange\"",
            "[ssa]\nrate_correction = true",
            "[ssa]\nspecies_levels = 2",
            "[ssa.species_levels]\nO2 = -1",
            "[ssa.species_levels]\nO2 = 1.5",
            "[ssa.species_levels]\nO2 = \"2\"",
            "[ssa.species_levels]\nO2 = 17",
        ] {
            let e = SolverConfig::from_toml_str(bad).unwrap_err().to_string();
            assert!(e.contains("solver") || e.contains("ssa"), "{bad}: {e}");
        }
    }

    #[test]
    fn unit_cell_refuses_an_explicit_mesh_correction() {
        let parse = |extra: &str| {
            SolverConfig::from_toml_str(&format!("[ssa]\nvolume = 1e-15\n{extra}"))
                .unwrap()
                .ssa
        };
        assert_eq!(parse("").require_volume().unwrap(), 1e-15);
        assert_eq!(
            parse("rate_correction = \"none\"")
                .require_volume()
                .unwrap(),
            1e-15
        );
        let e = parse("rate_correction = \"hellander\"")
            .require_volume()
            .unwrap_err();
        assert!(e.to_string().contains("rate_correction"), "{e}");
    }

    #[test]
    fn parses_background_run() {
        let parse = |extra: &str| {
            SolverConfig::from_toml_str(&format!(
                "[ssa]
{extra}"
            ))
            .map(|c| c.ssa)
        };
        assert_eq!(parse("").unwrap().background_run, None);
        assert_eq!(
            parse("background_run = false").unwrap().background_run,
            Some(false)
        );
        assert_eq!(
            parse("background_run = true").unwrap().background_run,
            Some(true)
        );
        for bad in ["background_run = \"yes\"", "background_run = 1"] {
            let e = parse(bad).unwrap_err();
            assert!(e.to_string().contains("background_run"), "{e}");
        }
        let e = parse(
            "volume = 1e-15
background_run = true",
        )
        .unwrap()
        .require_volume()
        .unwrap_err();
        assert!(e.to_string().contains("background_run"), "{e}");
    }

    #[test]
    fn parses_species_levels() {
        let c = SolverConfig::from_toml_str("[ssa.species_levels]\nO2 = 2\nH3O_plus = 0")
            .unwrap()
            .ssa;
        assert_eq!(
            c.species_levels,
            [("O2".to_string(), 2), ("H3O_plus".to_string(), 0)]
        );
        let species = ["e_aq", "O2", "OH_r"].map(String::from);
        let e = c.resolve_species_levels(&species).unwrap_err().to_string();
        assert!(e.contains("ssa.species_levels.H3O_plus"), "{e}");
        let c = SolverConfig::from_toml_str("[ssa.species_levels]\nO2 = 2")
            .unwrap()
            .ssa;
        assert_eq!(c.resolve_species_levels(&species).unwrap(), [0, 2, 0]);
        let e = SolverConfig::from_toml_str("[ssa.species_levels]\nO2 = -1")
            .unwrap_err()
            .to_string();
        assert!(e.contains("ssa.species_levels.O2"), "{e}");
    }

    #[test]
    fn species_levels_need_a_stochastic_voxel_run() {
        let c = SolverConfig::from_toml_str("[ssa]\nvolume = 1e-15\n[ssa.species_levels]\nO2 = 1")
            .unwrap()
            .ssa;
        let e = c.require_volume().unwrap_err().to_string();
        assert!(
            e.contains("species_levels") && e.contains("unit-cell"),
            "{e}"
        );
        let e = c
            .refuse_species_levels("a deterministic voxel run")
            .unwrap_err()
            .to_string();
        assert!(e.contains("deterministic"), "{e}");
        assert!(SsaConfig::default().refuse_species_levels("x").is_ok());
    }

    #[test]
    fn depth_covers_coarsening_and_species() {
        assert_eq!(hierarchy_depth(0, &[]), 0);
        assert_eq!(hierarchy_depth(2, &[0, 1]), 2);
        assert_eq!(hierarchy_depth(1, &[0, 3, 2]), 3);
    }

    #[test]
    fn suggestion_follows_the_square_root_of_the_diffusion_ratio() {
        // sqrt(16) = 4 = 2^2.
        assert_eq!(suggested_level(16e-9, 1e-9, 5), Some(2));
        // sqrt(15) < 4: one level.
        assert_eq!(suggested_level(15e-9, 1e-9, 5), Some(1));
        // Capped by the depth, never below 0.
        assert_eq!(suggested_level(1e-6, 1e-9, 2), Some(2));
        assert_eq!(suggested_level(1e-10, 1e-9, 2), Some(0));
        assert_eq!(suggested_level(0.0, 1e-9, 2), None);
        assert_eq!(suggested_level(1e-9, 0.0, 2), None);
    }
}
