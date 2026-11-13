//! Named beams with their LET (port of radiopyo `physics/beam_collection.py`).

use indexmap::IndexMap;

use crate::error::{RadkinError, Result};
use crate::physics::beam::Beam;

/// Beams by label, each with an optional LET [keV/µm], and the current one.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BeamCollection {
    beams: IndexMap<String, (Beam, Option<f64>)>,
    current: String,
}

fn number(v: &toml::Value) -> Option<f64> {
    match v {
        toml::Value::Float(f) => Some(*f),
        toml::Value::Integer(i) => Some(*i as f64),
        _ => None,
    }
}

impl BeamCollection {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.beams.len()
    }

    pub fn is_empty(&self) -> bool {
        self.beams.is_empty()
    }

    pub fn contains(&self, key: &str) -> bool {
        self.beams.contains_key(key)
    }

    /// Labels in insertion order.
    pub fn labels(&self) -> Vec<&str> {
        self.beams.keys().map(String::as_str).collect()
    }

    pub fn get(&self, key: &str) -> Option<&(Beam, Option<f64>)> {
        self.beams.get(key)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &(Beam, Option<f64>))> {
        self.beams.iter()
    }

    /// Current beam.
    pub fn current(&self) -> Result<&Beam> {
        self.beams
            .get(&self.current)
            .map(|(b, _)| b)
            .ok_or(RadkinError::NoBeamDefined)
    }

    pub fn current_name(&self) -> &str {
        &self.current
    }

    /// LET of the current beam.
    pub fn current_let(&self) -> Result<f64> {
        match self.beams.get(&self.current) {
            Some((_, Some(l))) => Ok(*l),
            _ => Err(RadkinError::Value(
                "LET not defined for current beam".into(),
            )),
        }
    }

    /// LET of the current beam, `None` if unset or no beam.
    pub fn current_let_opt(&self) -> Option<f64> {
        self.beams.get(&self.current).and_then(|(_, l)| *l)
    }

    pub fn remove(&mut self, key: &str) -> Result<()> {
        self.beams
            .shift_remove(key)
            .map(|_| ())
            .ok_or_else(|| RadkinError::Beam(format!("No beam definition named: {key}")))
    }

    /// Set the LET of `key` (current beam if `None`).
    pub fn update_let(&mut self, let_: Option<f64>, key: Option<&str>) -> Result<()> {
        let key = key.unwrap_or(&self.current).to_string();
        match self.beams.get_mut(&key) {
            Some(entry) => {
                entry.1 = let_;
                Ok(())
            }
            None => Err(RadkinError::Beam(format!(
                "No beam definition named: {key}"
            ))),
        }
    }

    /// Add (or replace, moved to the end) a beam.
    pub fn add_beam(&mut self, beam: Beam, key: &str, let_: Option<f64>, use_: bool) {
        let key = if key.is_empty() { "default" } else { key };
        self.beams.shift_remove(key);
        self.beams.insert(key.to_string(), (beam, let_));
        if use_ || self.current.is_empty() {
            self.current = key.to_string();
        }
    }

    pub fn use_beam(&mut self, key: &str) -> Result<()> {
        if !self.beams.contains_key(key) {
            return Err(RadkinError::Beam(format!(
                "No beam definition named: {key}"
            )));
        }
        self.current = key.to_string();
        Ok(())
    }

    /// radiopyo `add_beam_from_dict`: LET and label read from the table.
    pub fn add_beam_from_table(
        &mut self,
        beam_type: &str,
        config: &toml::Table,
        let_: Option<f64>,
    ) -> Result<()> {
        let let_ = match let_ {
            Some(l) => Some(l),
            None => config.get("LET").and_then(number),
        };
        let label = match config.get("label") {
            None => "default".to_string(),
            Some(toml::Value::String(s)) => s.clone(),
            Some(v) => v.to_string(),
        };
        let beam = Beam::from_table(beam_type, config)?;
        self.add_beam(beam, &label, let_, false);
        Ok(())
    }

    /// Beams of a config `[beam]` section (`[[beam.<type>]]` tables).
    pub fn from_config(section: &toml::Table) -> Result<Self> {
        let mut out = Self::new();
        for (beam_type, value) in section {
            let tables: Vec<&toml::Table> = match value {
                toml::Value::Array(a) => a
                    .iter()
                    .map(|v| {
                        v.as_table().ok_or_else(|| {
                            RadkinError::Beam(format!("beam.{beam_type}: not a table"))
                        })
                    })
                    .collect::<Result<_>>()?,
                toml::Value::Table(t) => vec![t],
                _ => return Err(RadkinError::Beam(format!("beam.{beam_type}: not a table"))),
            };
            for t in tables {
                let mut config = t.clone();
                let let_ = config.remove("LET").as_ref().and_then(number);
                out.add_beam_from_table(beam_type, &config, let_)?;
            }
        }
        Ok(out)
    }

    /// Back to a `[beam]` section (radiopyo `as_dict`).
    pub fn as_config(&self) -> toml::Table {
        let mut out = toml::Table::new();
        for (name, (beam, let_)) in &self.beams {
            let mut x = beam.as_table();
            if name != "default" {
                x.insert("label".into(), toml::Value::String(name.clone()));
            }
            if let Some(l) = let_ {
                x.insert("LET".into(), toml::Value::Float(*l));
            }
            let entry = out
                .entry(beam.beam_type())
                .or_insert_with(|| toml::Value::Array(Vec::new()));
            if let toml::Value::Array(a) = entry {
                a.push(toml::Value::Table(x));
            }
        }
        out
    }

    /// radiopyo `BeamCollection.merge`: same label and type -> parameters
    /// updated, same label and another type -> replaced.
    pub fn merge(&self, other: &Self) -> Result<Self> {
        let mut s: IndexMap<String, (toml::Table, String)> = IndexMap::new();
        for (name, (beam, let_)) in &self.beams {
            let mut d = beam.as_table();
            if let Some(l) = let_ {
                d.insert("LET".into(), toml::Value::Float(*l));
            }
            d.insert("label".into(), toml::Value::String(name.clone()));
            s.insert(name.clone(), (d, beam.beam_type().to_string()));
        }
        for (beam_type, beams) in other.as_config() {
            for beam in beams
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(|b| b.as_table())
            {
                let name = match beam.get("label") {
                    Some(toml::Value::String(l)) => l.clone(),
                    _ => "default".to_string(),
                };
                match s.get_mut(&name) {
                    None => {
                        s.insert(name, (beam.clone(), beam_type.clone()));
                    }
                    Some((d, t)) if *t == beam_type => {
                        for (k, v) in beam {
                            d.insert(k.clone(), v.clone());
                        }
                    }
                    Some(entry) => *entry = (beam.clone(), beam_type.clone()),
                }
            }
        }
        let mut out = Self::new();
        for (_, (beam, beam_type)) in s {
            out.add_beam_from_table(&beam_type, &beam, None)?;
        }
        Ok(out)
    }
}
