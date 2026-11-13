//! Reaction matrices built once per environment (port of radiopyo
//! `SimEnv.build_all_matrices` and `utils/sim_types.py`).
//!
//! Index conventions: species `0..n_dyn` are dynamic, `n_dyn..n_species`
//! fixed; the padding index `n_species` stands for a factor 1.

use crate::chemistry::env::SimEnv;
use crate::chemistry::parser::Operation;
use crate::error::Result;

/// Mass-action reactions (k-reactions, acid/base and equilibrium pairs).
#[derive(Debug, Clone, PartialEq)]
pub struct KMatrices {
    pub k_values: Vec<f64>,
    /// Reactant species of each reaction, one slot per unit of stoichiometry.
    pub reactant_index: Vec<Vec<usize>>,
    /// Non-zero stoichiometry `(dyn species, coef)` of each reaction.
    pub stoichiometry: Vec<Vec<(usize, f64)>>,
    /// "Negative" product species of each reaction (depletion switch).
    pub negative_index: Vec<Vec<usize>>,
}

/// Radiolytic production.
#[derive(Debug, Clone, PartialEq)]
pub struct GVector {
    /// Yields [mol/L/Gy].
    pub yields: Vec<f64>,
    /// Produced dynamic species (`None` if the species is fixed).
    pub species: Vec<Option<usize>>,
}

/// Michaelis–Menten reactions.
#[derive(Debug, Clone, PartialEq)]
pub struct EzMatrices {
    pub k_values: Vec<f64>,
    pub k_michaelis: Vec<f64>,
    pub substrate_index: Vec<usize>,
    pub enzyme_index: Vec<usize>,
    pub stoichiometry: Vec<Vec<(usize, f64)>>,
    pub negative_index: Vec<Vec<usize>>,
}

/// One mass-action reaction seen from a dynamic species `j` (overshoot path).
#[derive(Debug, Clone, PartialEq)]
pub struct KOvershootEntry {
    pub k: f64,
    pub slots: Vec<usize>,
    /// `k` if `j` is exactly one reactant slot, else 0.
    pub k_dj: f64,
    /// `slots` with `j` replaced by the padding index.
    pub slots_dj: Vec<usize>,
    /// Stoichiometry of `j` in the reaction.
    pub weight: f64,
    /// `j` is a negative product of the reaction.
    pub flag: bool,
    /// Negative products other than `j`.
    pub negative_without_j: Vec<usize>,
}

/// One enzymatic reaction seen from a dynamic species `j`.
#[derive(Debug, Clone, PartialEq)]
pub struct EzOvershootEntry {
    pub k: f64,
    pub km: f64,
    pub substrate: usize,
    pub enzyme: usize,
    pub is_enzyme: bool,
    pub is_substrate: bool,
    pub weight: f64,
    pub flag: bool,
    pub negative_without_j: Vec<usize>,
}

/// Per dynamic species tables of the overshoot relaxation (ADR 0007).
#[derive(Debug, Clone, PartialEq)]
pub struct OvershootTables {
    pub k: Vec<Vec<KOvershootEntry>>,
    pub ez: Vec<Vec<EzOvershootEntry>>,
}

/// Every matrix of the ODE system.
#[derive(Debug, Clone, PartialEq)]
pub struct SimMatrices {
    pub k: KMatrices,
    pub radiolytic: GVector,
    pub enzymatic: EzMatrices,
    /// Fixed species concentrations [mol/L], species `n_dyn..`.
    pub cst_sp: Vec<f64>,
    pub overshoot: OvershootTables,
    pub n_species: usize,
    pub n_dyn: usize,
}

fn stoichiometry_entry(col: &mut Vec<(usize, f64)>, species: usize, delta: f64) {
    match col.iter_mut().find(|(s, _)| *s == species) {
        Some(e) => e.1 += delta,
        None => col.push((species, delta)),
    }
}

fn weight_of(col: &[(usize, f64)], j: usize) -> f64 {
    col.iter()
        .find(|(s, _)| *s == j)
        .map(|(_, c)| *c)
        .unwrap_or(0.0)
}

impl SimMatrices {
    /// Build every matrix at the environment temperature.
    pub fn build(env: &SimEnv) -> Result<Self> {
        let sp = &env.species;
        let n_species = sp.len();
        let n_dyn = sp.n_dyn();
        let idx = |l: &str| sp.index_of(l);
        let mut k = KMatrices {
            k_values: Vec::new(),
            reactant_index: Vec::new(),
            stoichiometry: Vec::new(),
            negative_index: Vec::new(),
        };
        for r in env.reactions.iter_all_k() {
            k.k_values.push(r.k_value(env.temperature)?);
            let mut slots = Vec::new();
            let mut col = Vec::new();
            let mut negative = Vec::new();
            for s in &r.reactants {
                let i = idx(&s.label)?;
                if !sp.is_cst(&s.label) {
                    stoichiometry_entry(&mut col, i, -(s.stoi as f64));
                }
                slots.extend(std::iter::repeat_n(i, s.stoi as usize));
            }
            for s in &r.products {
                if sp.is_cst(&s.label) {
                    continue;
                }
                let i = idx(&s.label)?;
                if s.ope == Operation::Add {
                    stoichiometry_entry(&mut col, i, s.stoi as f64);
                } else {
                    stoichiometry_entry(&mut col, i, -(s.stoi as f64));
                    negative.push(i);
                }
            }
            col.retain(|(_, c)| *c != 0.0);
            negative.sort_unstable();
            k.reactant_index.push(slots);
            k.stoichiometry.push(col);
            k.negative_index.push(negative);
        }
        let mut radiolytic = GVector {
            yields: Vec::new(),
            species: Vec::new(),
        };
        for r in &env.reactions.radiolytic {
            radiolytic.yields.push(r.kr);
            radiolytic.species.push(if sp.is_cst(&r.species) {
                None
            } else {
                Some(idx(&r.species)?)
            });
        }
        let mut ez = EzMatrices {
            k_values: Vec::new(),
            k_michaelis: Vec::new(),
            substrate_index: Vec::new(),
            enzyme_index: Vec::new(),
            stoichiometry: Vec::new(),
            negative_index: Vec::new(),
        };
        for r in &env.reactions.enzymatic {
            ez.substrate_index.push(idx(r.substrate())?);
            ez.enzyme_index.push(idx(&r.enzyme)?);
            ez.k_michaelis.push(r.k_micha);
            ez.k_values.push(r.k_value);
            let mut col = Vec::new();
            let mut negative = Vec::new();
            for s in &r.reactants {
                if !sp.is_cst(&s.label) {
                    stoichiometry_entry(&mut col, idx(&s.label)?, -(s.stoi as f64));
                }
            }
            for s in &r.products {
                if sp.is_cst(&s.label) {
                    continue;
                }
                let i = idx(&s.label)?;
                if s.ope == Operation::Add {
                    stoichiometry_entry(&mut col, i, s.stoi as f64);
                } else {
                    stoichiometry_entry(&mut col, i, -(s.stoi as f64));
                    negative.push(i);
                }
            }
            col.retain(|(_, c)| *c != 0.0);
            negative.sort_unstable();
            ez.stoichiometry.push(col);
            ez.negative_index.push(negative);
        }
        let mut cst_sp = vec![0.0; sp.n_cst()];
        for s in sp.cst_sp.values() {
            cst_sp[s.index - n_dyn] = s.initial_cc(env.temperature);
        }
        let overshoot = OvershootTables::build(&k, &ez, n_dyn, n_species);
        Ok(Self {
            k,
            radiolytic,
            enzymatic: ez,
            cst_sp,
            overshoot,
            n_species,
            n_dyn,
        })
    }
}

impl OvershootTables {
    fn build(k: &KMatrices, ez: &EzMatrices, n_dyn: usize, n_species: usize) -> Self {
        let mut kt = Vec::with_capacity(n_dyn);
        let mut et = Vec::with_capacity(n_dyn);
        for j in 0..n_dyn {
            let mut entries = Vec::new();
            for (r, col) in k.stoichiometry.iter().enumerate() {
                let weight = weight_of(col, j);
                if weight == 0.0 {
                    continue;
                }
                let slots = k.reactant_index[r].clone();
                let once = slots.iter().filter(|&&s| s == j).count() == 1;
                entries.push(KOvershootEntry {
                    k: k.k_values[r],
                    k_dj: if once { k.k_values[r] } else { 0.0 },
                    slots_dj: slots
                        .iter()
                        .map(|&s| if s == j { n_species } else { s })
                        .collect(),
                    slots,
                    weight,
                    flag: k.negative_index[r].contains(&j),
                    negative_without_j: k.negative_index[r]
                        .iter()
                        .copied()
                        .filter(|&s| s != j)
                        .collect(),
                });
            }
            kt.push(entries);
            let mut entries = Vec::new();
            for (r, col) in ez.stoichiometry.iter().enumerate() {
                let weight = weight_of(col, j);
                if weight == 0.0 {
                    continue;
                }
                entries.push(EzOvershootEntry {
                    k: ez.k_values[r],
                    km: ez.k_michaelis[r],
                    substrate: ez.substrate_index[r],
                    enzyme: ez.enzyme_index[r],
                    is_enzyme: ez.enzyme_index[r] == j,
                    is_substrate: ez.substrate_index[r] == j,
                    weight,
                    flag: ez.negative_index[r].contains(&j),
                    negative_without_j: ez.negative_index[r]
                        .iter()
                        .copied()
                        .filter(|&s| s != j)
                        .collect(),
                });
            }
            et.push(entries);
        }
        Self { k: kt, ez: et }
    }
}
