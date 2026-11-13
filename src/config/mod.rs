//! TOML configuration: loading, `[includes]` resolution and merge rules
//! (port of radiopyo `parser/toml.py`, `parser/file_parser.py`,
//! `parser/utils.py`).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;

use crate::chemistry::parser::{ReactionKey, reaction_key};
use crate::consts::{DEFAULT_TEMPERATURE, data_dir};
use crate::error::{RadkinError, Result};
use crate::physics::BeamCollection;

/// Prefixes of include entries that refer to the bundled data folder.
pub const PACKAGE_PREFIXES: [&str; 2] = ["radiopyo/", "radkin/"];

/// Reaction list sections of `[reactions]`, in radiopyo merge order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Section {
    AcidBase,
    Equilibrium,
    KReaction,
    Constant,
    Arrhenius,
    Enzymatic,
}

impl Section {
    pub const ALL: [Section; 6] = [
        Section::AcidBase,
        Section::Equilibrium,
        Section::KReaction,
        Section::Constant,
        Section::Arrhenius,
        Section::Enzymatic,
    ];

    /// TOML name of the section.
    pub fn name(&self) -> &'static str {
        match self {
            Self::AcidBase => "acid_base",
            Self::Equilibrium => "equilibrium",
            Self::KReaction => "k_reaction",
            Self::Constant => "constant",
            Self::Arrhenius => "arrhenius",
            Self::Enzymatic => "enzymatic",
        }
    }
}

/// A configuration (radiopyo `ConfigDict`).
#[derive(Debug, Clone, PartialEq)]
pub struct Config {
    /// `ext -> entries` of `[includes]`, anchored to the file folder.
    pub includes: IndexMap<String, Vec<String>>,
    pub parameters: toml::Table,
    /// `[concentrations.fixed]` [mol/L] (number or temperature polynomial).
    pub fixed: toml::Table,
    /// `[concentrations.initial]` [mol/L].
    pub initial: toml::Table,
    /// `[reactions.radiolytic]` G-values [#/100 eV].
    pub radiolytic: toml::Table,
    /// Every reaction list section (always present).
    pub reactions: IndexMap<Section, Vec<toml::Table>>,
    /// `[beam]` section (`type -> [tables]`).
    pub beam: toml::Table,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            includes: IndexMap::new(),
            parameters: toml::Table::new(),
            fixed: toml::Table::new(),
            initial: toml::Table::new(),
            radiolytic: toml::Table::new(),
            reactions: Section::ALL.iter().map(|s| (*s, Vec::new())).collect(),
            beam: toml::Table::new(),
        }
    }
}

fn as_number(v: &toml::Value) -> Option<f64> {
    match v {
        toml::Value::Float(f) => Some(*f),
        toml::Value::Integer(i) => Some(*i as f64),
        _ => None,
    }
}

fn sub_table(t: &toml::Table, key: &str, what: &str) -> Result<toml::Table> {
    match t.get(key) {
        None => Ok(toml::Table::new()),
        Some(toml::Value::Table(x)) => Ok(x.clone()),
        Some(_) => Err(RadkinError::Config(format!("'{what}' must be a table"))),
    }
}

fn table_list(v: Option<&toml::Value>, what: &str) -> Result<Vec<toml::Table>> {
    match v {
        None => Ok(Vec::new()),
        Some(toml::Value::Table(t)) => Ok(if t.is_empty() {
            Vec::new()
        } else {
            vec![t.clone()]
        }),
        Some(toml::Value::Array(a)) => a
            .iter()
            .map(|x| {
                x.as_table()
                    .cloned()
                    .ok_or_else(|| RadkinError::Config(format!("'{what}' entries must be tables")))
            })
            .collect(),
        Some(_) => Err(RadkinError::Config(format!(
            "'{what}' must be a list of tables"
        ))),
    }
}

/// Path of an include entry (radiopyo `include_path`).
pub fn include_path(entry: &str) -> PathBuf {
    let entry = entry.trim();
    for prefix in PACKAGE_PREFIXES {
        if let Some(rest) = entry.strip_prefix(prefix) {
            return data_dir().join(rest);
        }
    }
    std::path::absolute(entry).unwrap_or_else(|_| PathBuf::from(entry))
}

fn anchored(entry: &str, base: &Path) -> String {
    let entry = entry.trim();
    if PACKAGE_PREFIXES.iter().any(|p| entry.starts_with(p)) || Path::new(entry).is_absolute() {
        return entry.to_string();
    }
    let candidate = base.join(entry);
    if candidate.exists() {
        candidate.to_string_lossy().into_owned()
    } else {
        entry.to_string()
    }
}

impl Config {
    /// Parse one TOML file, without resolving its includes.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
        Self::from_toml_str(&text, path.parent())
            .map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
    }

    /// Parse TOML text; relative includes are anchored to `base`.
    pub fn from_toml_str(text: &str, base: Option<&Path>) -> Result<Self> {
        let root: toml::Table =
            toml::from_str(text).map_err(|e| RadkinError::Config(e.to_string()))?;
        Self::from_table(root, base)
    }

    /// Configuration from an already parsed TOML table (e.g. a radiopyo
    /// `ConfigDict` built in Python); relative includes are anchored to `base`.
    pub fn from_table(root: toml::Table, base: Option<&Path>) -> Result<Self> {
        let mut out = Self::default();
        let base = std::path::absolute(base.unwrap_or(Path::new("."))).unwrap_or_default();
        for (ext, entries) in sub_table(&root, "includes", "includes")? {
            let list: Vec<String> = match entries {
                toml::Value::String(s) => vec![s],
                toml::Value::Array(a) => a
                    .into_iter()
                    .map(|v| match v {
                        toml::Value::String(s) => Ok(s),
                        _ => Err(RadkinError::Config(
                            "include entries must be strings".into(),
                        )),
                    })
                    .collect::<Result<_>>()?,
                _ => {
                    return Err(RadkinError::Config(
                        "include entries must be strings".into(),
                    ));
                }
            };
            out.includes
                .insert(ext, list.iter().map(|e| anchored(e, &base)).collect());
        }
        out.parameters = sub_table(&root, "parameters", "parameters")?;
        out.beam = sub_table(&root, "beam", "beam")?;
        let conc = sub_table(&root, "concentrations", "concentrations")?;
        for key in conc.keys() {
            if key != "fixed" && key != "initial" {
                log::warn!(
                    "CONFIG::Found unknown, extra, key(s) in 'concentrations' section: {key}. These will not be used."
                );
            }
        }
        out.fixed = sub_table(&conc, "fixed", "concentrations.fixed")?;
        out.initial = sub_table(&conc, "initial", "concentrations.initial")?;
        let reactions = sub_table(&root, "reactions", "reactions")?;
        out.radiolytic = sub_table(&reactions, "radiolytic", "reactions.radiolytic")?;
        for s in Section::ALL {
            out.reactions
                .insert(s, table_list(reactions.get(s.name()), s.name())?);
        }
        Ok(out)
    }

    /// Load a file and resolve its includes recursively.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        resolve_includes(Self::from_toml_file(path)?, &[])
    }

    fn param(&self, key: &str) -> Option<&toml::Value> {
        self.parameters.get(key)
    }

    /// Initial pH (default 7).
    pub fn ph(&self) -> f64 {
        self.param("pH").and_then(as_number).unwrap_or(7.0)
    }

    pub fn force_cst_ph(&self) -> bool {
        self.param("force_cst_pH")
            .and_then(|v| v.as_bool())
            .unwrap_or(false)
    }

    /// Temperature [degC] (default 25).
    pub fn temperature(&self) -> f64 {
        self.param("temperature")
            .and_then(as_number)
            .unwrap_or(DEFAULT_TEMPERATURE)
    }

    pub fn variable_gvalues(&self) -> bool {
        self.param("variable_GValues")
            .and_then(|v| v.as_bool())
            .unwrap_or(false)
    }

    /// Explicit `(dynamical, constant)` species order, if given.
    pub fn species_order(&self) -> Option<(Vec<String>, Vec<String>)> {
        let t = self.param("species_order")?.as_table()?;
        let list = |k: &str| -> Option<Vec<String>> {
            t.get(k)?
                .as_array()?
                .iter()
                .map(|v| v.as_str().map(String::from))
                .collect()
        };
        let dynamical = list("dynamical")?;
        if dynamical.is_empty() {
            return None;
        }
        Some((dynamical, list("constant").unwrap_or_default()))
    }
}

fn update(into: &mut toml::Table, from: &toml::Table) {
    for (k, v) in from {
        into.insert(k.clone(), v.clone());
    }
}

const EQUILIBRIUM_CONSTANTS: [&str; 3] = ["k_forward", "k_backward", "k_equilibrium"];

fn reaction_string(r: &toml::Table, section: Section) -> Result<&str> {
    r.get("reaction").and_then(|v| v.as_str()).ok_or_else(|| {
        RadkinError::Config(format!(
            "{}: entry without a 'reaction' string",
            section.name()
        ))
    })
}

fn merge_reaction_list(
    section: Section,
    base: &[toml::Table],
    over: &[toml::Table],
) -> Result<Vec<toml::Table>> {
    let mut position: HashMap<ReactionKey, usize> = HashMap::new();
    for (i, r) in base.iter().enumerate() {
        position.insert(reaction_key(reaction_string(r, section)?)?, i);
    }
    let mut merged = base.to_vec();
    for entry in over {
        let key = reaction_key(reaction_string(entry, section)?)?;
        let Some(&idx) = position.get(&key) else {
            merged.push(entry.clone());
            continue;
        };
        log::debug!(
            "CONFIG::{}: reaction redefined by '{}'",
            section.name(),
            reaction_string(entry, section)?
        );
        if section == Section::Equilibrium {
            let given = EQUILIBRIUM_CONSTANTS
                .iter()
                .filter(|c| entry.contains_key(**c))
                .count();
            if given >= 2 {
                for c in EQUILIBRIUM_CONSTANTS {
                    merged[idx].remove(c);
                }
            }
        }
        update(&mut merged[idx], entry);
    }
    Ok(merged)
}

/// Merge `other` over `root` (radiopyo `combine_config`).
pub fn combine(root: &Config, other: &Config) -> Result<Config> {
    let mut out = root.clone();
    out.includes.clear();
    update(&mut out.parameters, &other.parameters);
    update(&mut out.fixed, &other.fixed);
    update(&mut out.initial, &other.initial);
    let mut root_beams = BeamCollection::from_config(&out.beam)?;
    let other_beams = BeamCollection::from_config(&other.beam)?;
    if root_beams.contains("default") && !other_beams.is_empty() {
        root_beams.remove("default")?;
    }
    out.beam = root_beams.merge(&other_beams)?.as_config();
    update(&mut out.radiolytic, &other.radiolytic);
    for s in Section::ALL {
        let merged = merge_reaction_list(s, &out.reactions[&s], &other.reactions[&s])?;
        out.reactions.insert(s, merged);
    }
    Ok(out)
}

/// Resolve the includes of `config` recursively (radiopyo `resolve_includes`).
pub fn resolve_includes(config: Config, chain: &[PathBuf]) -> Result<Config> {
    if config.includes.is_empty() {
        return Ok(config);
    }
    let mut merged = Config::default();
    for (ext, files) in &config.includes {
        if ext.trim_start_matches('.') != "toml" {
            log::warn!("Unable to include .{ext} file. No parser available");
            continue;
        }
        for f in files {
            let path = include_path(f);
            if !path.exists() {
                return Err(RadkinError::Config(format!(
                    "include not found: {}",
                    path.display()
                )));
            }
            let resolved = path.canonicalize().map_err(|e| RadkinError::io(&path, e))?;
            if chain.contains(&resolved) {
                let cycle: Vec<String> = chain
                    .iter()
                    .chain([&resolved])
                    .map(|p| p.display().to_string())
                    .collect();
                return Err(RadkinError::Config(format!(
                    "Circular include: {}",
                    cycle.join(" -> ")
                )));
            }
            log::info!("Resolving include: {}", path.display());
            let mut sub_chain = chain.to_vec();
            sub_chain.push(resolved);
            let included = resolve_includes(Config::from_toml_file(&path)?, &sub_chain)?;
            merged = combine(&merged, &included)?;
        }
    }
    combine(&merged, &config)
}
