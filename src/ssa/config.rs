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
//! ```
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
}

impl Default for SsaConfig {
    fn default() -> Self {
        Self {
            volume: None,
            trajectories: DEFAULT_TRAJECTORIES,
            seed: 0,
            tau_leap: false,
            tau_epsilon: DEFAULT_TAU_EPSILON,
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
                other => {
                    return Err(RadkinError::Config(format!(
                        "unknown key ssa.{other} (known: volume, trajectories, seed, tau_leap, tau_epsilon)"
                    )));
                }
            }
        }
        Ok(out)
    }

    /// The unit-cell volume [L], an error if `volume` is absent.
    pub fn require_volume(&self) -> Result<f64> {
        self.volume.ok_or_else(|| {
            RadkinError::Config(
                "ssa.volume (litres) is required by a unit-cell run with solver = \"ssa\"".into(),
            )
        })
    }
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
tau_epsilon = 0.01",
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
        ] {
            let e = SolverConfig::from_toml_str(bad).unwrap_err().to_string();
            assert!(e.contains("solver") || e.contains("ssa"), "{bad}: {e}");
        }
    }
}
