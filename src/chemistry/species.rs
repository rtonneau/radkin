//! Simulation species (port of radiopyo `species.py` / `species_collection.py`).

use indexmap::IndexMap;

use crate::error::{RadkinError, Result};
use crate::poly::Polynomial;

/// A species of the simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct SimSpecies {
    pub label: String,
    /// Index in the full state vector (dynamic species first).
    pub index: usize,
    /// Initial (or fixed) concentration [mol/L] as a polynomial of the
    /// temperature [degC]; `None` = 0.
    pub cc_init: Option<Polynomial>,
}

impl SimSpecies {
    /// Initial concentration [mol/L] at `t_celsius` (negative values give 0).
    pub fn initial_cc(&self, t_celsius: f64) -> f64 {
        match &self.cc_init {
            None => 0.0,
            Some(p) => {
                let cc = p.eval(t_celsius);
                if cc > 0.0 { cc } else { 0.0 }
            }
        }
    }
}

/// Dynamic (solved) and constant (fixed) species, in order of appearance.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ListSpecies {
    pub dyn_sp: IndexMap<String, SimSpecies>,
    pub cst_sp: IndexMap<String, SimSpecies>,
    next_index: usize,
}

impl ListSpecies {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn contains(&self, label: &str) -> bool {
        self.dyn_sp.contains_key(label) || self.cst_sp.contains_key(label)
    }

    pub fn get(&self, label: &str) -> Option<&SimSpecies> {
        self.dyn_sp.get(label).or_else(|| self.cst_sp.get(label))
    }

    pub fn get_mut(&mut self, label: &str) -> Option<&mut SimSpecies> {
        match self.dyn_sp.get_mut(label) {
            Some(s) => Some(s),
            None => self.cst_sp.get_mut(label),
        }
    }

    /// Index of `label` in the full state vector.
    pub fn index_of(&self, label: &str) -> Result<usize> {
        self.get(label)
            .map(|s| s.index)
            .ok_or_else(|| RadkinError::Value(format!("unknown species '{label}'")))
    }

    pub fn is_cst(&self, label: &str) -> bool {
        self.cst_sp.contains_key(label)
    }

    pub fn is_dyn(&self, label: &str) -> bool {
        self.dyn_sp.contains_key(label)
    }

    pub fn n_dyn(&self) -> usize {
        self.dyn_sp.len()
    }

    pub fn n_cst(&self) -> usize {
        self.cst_sp.len()
    }

    pub fn len(&self) -> usize {
        self.n_dyn() + self.n_cst()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn add_dyn(&mut self, label: &str, cc_init: Option<Polynomial>) -> Result<()> {
        if self.contains(label) {
            return Err(RadkinError::Config(format!("{label} is already recorded")));
        }
        let sp = SimSpecies {
            label: label.to_string(),
            index: self.next_index,
            cc_init,
        };
        self.dyn_sp.insert(label.to_string(), sp);
        self.next_index += 1;
        Ok(())
    }

    pub fn add_cst(&mut self, label: &str, cc: Polynomial) -> Result<()> {
        if self.contains(label) {
            return Err(RadkinError::Config(format!("{label} is already recorded")));
        }
        let sp = SimSpecies {
            label: label.to_string(),
            index: self.next_index,
            cc_init: Some(cc),
        };
        self.cst_sp.insert(label.to_string(), sp);
        self.next_index += 1;
        Ok(())
    }

    /// Existing species (dynamic or constant) or a new dynamic one.
    pub fn get_or_create_dyn(&mut self, label: &str) {
        if !self.contains(label) {
            self.add_dyn(label, None).expect("checked");
        }
    }

    pub fn remove(&mut self, label: &str) {
        if self.dyn_sp.shift_remove(label).is_none() {
            self.cst_sp.shift_remove(label);
        }
    }

    /// radiopyo `ListSpecies.sort`: `starter` first, the other dynamic
    /// species in order of appearance, then the constant ones.
    pub fn sort(&mut self, starter: &[String]) -> Result<()> {
        for (idx, label) in starter.iter().enumerate() {
            if self.is_cst(label) {
                return Err(RadkinError::Config(format!(
                    "{label} is defined as constant. Not able to use it in an A/B equilibrium"
                )));
            }
            self.get_mut(label).expect("acid/base species exist").index = idx;
        }
        let others: Vec<String> = self
            .dyn_sp
            .keys()
            .filter(|l| !starter.contains(l))
            .cloned()
            .collect();
        for (idx, label) in others.iter().enumerate() {
            self.dyn_sp[label.as_str()].index = idx + starter.len();
        }
        let n_dyn = self.n_dyn();
        for (idx, sp) in self.cst_sp.values_mut().enumerate() {
            sp.index = idx + n_dyn;
        }
        self.next_index = self.len();
        Ok(())
    }

    /// radiopyo `ListSpecies.set_order`.
    pub fn set_order(&mut self, dyn_order: &[String], cst_order: &[String]) -> Result<()> {
        if dyn_order.len() != self.n_dyn() || cst_order.len() != self.n_cst() {
            return Err(RadkinError::Config(format!(
                "Unable to sort species, different numbers: {} dyn / {} cst read, order lists of {} / {}",
                self.n_dyn(),
                self.n_cst(),
                dyn_order.len(),
                cst_order.len()
            )));
        }
        for (idx, label) in dyn_order.iter().chain(cst_order).enumerate() {
            self.get_mut(label)
                .ok_or_else(|| {
                    RadkinError::Config(format!("species_order: unknown species '{label}'"))
                })?
                .index = idx;
        }
        self.next_index = self.len();
        Ok(())
    }

    /// Dynamic species sorted by index.
    pub fn sorted_dyn(&self) -> Vec<&SimSpecies> {
        let mut v: Vec<_> = self.dyn_sp.values().collect();
        v.sort_by_key(|s| s.index);
        v
    }

    /// Constant species sorted by index.
    pub fn sorted_cst(&self) -> Vec<&SimSpecies> {
        let mut v: Vec<_> = self.cst_sp.values().collect();
        v.sort_by_key(|s| s.index);
        v
    }

    pub fn sorted_dyn_labels(&self) -> Vec<String> {
        self.sorted_dyn()
            .into_iter()
            .map(|s| s.label.clone())
            .collect()
    }

    pub fn sorted_cst_labels(&self) -> Vec<String> {
        self.sorted_cst()
            .into_iter()
            .map(|s| s.label.clone())
            .collect()
    }
}
