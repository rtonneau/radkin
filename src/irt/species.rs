//! Species properties of the IRT stage: van der Waals radius, charge and
//! diffusion coefficient.
//!
//! The built-in table holds the water-radiolysis species under radkin's
//! labels, with the values Geant4-DNA 11.4.3 (chemistry option3) uses at run
//! time (`tests/data/irt/g4_chemistry_dump.txt`). A network can override or
//! add an entry.

use indexmap::IndexMap;

/// Properties of one species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpeciesProps {
    /// Van der Waals radius [m].
    pub radius: f64,
    /// Charge, in elementary charges.
    pub charge: i32,
    /// Diffusion coefficient [m^2/s].
    pub diffusion: f64,
}

const NM: f64 = 1e-9;

/// `(label, radius [nm], charge, D [m^2/s])` of the water-radiolysis species
/// (Geant4-DNA 11.4.3 option3).
pub const WATER_SPECIES: [(&str, f64, i32, f64); 15] = [
    ("H_r", 0.19, 0, 7.0e-9),
    ("H2", 0.14, 0, 4.8e-9),
    ("H2O2", 0.21, 0, 2.3e-9),
    ("H3O_plus", 0.25, 1, 9.46e-9),
    ("HO2_minus", 0.25, -1, 1.4e-9),
    ("HO2", 0.21, 0, 2.3e-9),
    ("O_atom", 0.2, 0, 2.0e-9),
    ("O2", 0.17, 0, 2.4e-9),
    ("O2_minus", 0.22, -1, 1.75e-9),
    ("O3", 0.2, 0, 2.0e-9),
    ("O3_minus", 0.2, -1, 2.0e-9),
    ("OH_minus", 0.33, -1, 5.3e-9),
    ("O_minus", 0.25, -1, 2.0e-9),
    ("e_aq", 0.5, -1, 4.9e-9),
    ("OH_r", 0.22, 0, 2.2e-9),
];

/// Species properties by label.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpeciesTable {
    entries: IndexMap<String, SpeciesProps>,
}

impl SpeciesTable {
    /// Table without any entry.
    pub fn empty() -> Self {
        Self::default()
    }

    /// The built-in table of the water-radiolysis species.
    pub fn water() -> Self {
        let entries = WATER_SPECIES
            .iter()
            .map(|&(label, radius, charge, diffusion)| {
                (
                    label.to_string(),
                    SpeciesProps {
                        radius: radius * NM,
                        charge,
                        diffusion,
                    },
                )
            })
            .collect();
        Self { entries }
    }

    /// Add or replace the entry of `label`.
    pub fn set(&mut self, label: &str, props: SpeciesProps) {
        self.entries.insert(label.to_string(), props);
    }

    pub fn get(&self, label: &str) -> Option<&SpeciesProps> {
        self.entries.get(label)
    }

    pub fn labels(&self) -> impl Iterator<Item = &str> {
        self.entries.keys().map(|k| k.as_str())
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_table_and_override() {
        let mut t = SpeciesTable::water();
        assert_eq!(t.len(), 15);
        let e = *t.get("e_aq").unwrap();
        assert_eq!((e.radius, e.charge, e.diffusion), (0.5e-9, -1, 4.9e-9));
        assert!(t.get("Cl_minus").is_none());
        t.set(
            "e_aq",
            SpeciesProps {
                radius: 0.3e-9,
                ..e
            },
        );
        t.set(
            "Cl_minus",
            SpeciesProps {
                radius: 0.18e-9,
                charge: -1,
                diffusion: 2.0e-9,
            },
        );
        assert_eq!(t.get("e_aq").unwrap().radius, 0.3e-9);
        assert_eq!(t.len(), 16);
    }
}
