//! Reaction model of the IRT stage, built from a network.
//!
//! Every mass-action reaction of [`SimMatrices`] is sorted by its dynamic
//! reactants:
//!
//! - two dynamic reactants: a **pair reaction** between two particles, with a
//!   type, its radii and its rates. If one reactant is also a background
//!   species (non-zero initial concentration), the other one gets a
//!   **molecule reaction** with that background as well;
//! - one dynamic reactant (first order, or with fixed species such as water):
//!   a **molecule reaction** with an exponential reaction time;
//! - no dynamic reactant: skipped (the background keeps itself at
//!   equilibrium outside the IRT stage).
//!
//! The pair formulas are those of Geant4-DNA 11.4.3
//! (`G4DNAMolecularReactionTable.cc`). `k_obs` is the rate constant of the
//! reaction events, `events = k_obs [A][B]` and `events = k_obs [A]^2` for
//! identical reactants: this is the value radkin stores for every reaction
//! (a scalar `k_value` of `2 A -> ...` is divided by 2, the n·k convention),
//! and it is Geant4's observed reaction rate.

use std::f64::consts::PI;

use crate::chemistry::env::SimEnv;
use crate::chemistry::parser::Operation;
use crate::chemistry::reactions::IrtType;
use crate::consts::AVOGADRO;
use crate::error::{RadkinError, Result};
use crate::irt::species::SpeciesTable;
use crate::matrices::SimMatrices;

/// Elementary charge [C] (CLHEP `e_SI`).
const E_SI: f64 = 1.602176634e-19;
/// Boltzmann constant [J/K], from CLHEP's `k_Boltzmann = 8.617333e-11 MeV/K`.
const K_BOLTZMANN: f64 = 8.617333e-5 * E_SI;
/// Vacuum permittivity [F/m], from CLHEP's `mu0 = 4 pi 1e-7 H/m`.
const EPSILON0: f64 = 1.0 / (4.0 * PI * 1e-7 * 299_792_458.0 * 299_792_458.0);
/// Temperature of the Onsager radius [K] (Geant4).
pub const ONSAGER_TEMPERATURE: f64 = 293.15;
/// Relative permittivity of water of the Onsager radius (Geant4).
pub const WATER_PERMITTIVITY: f64 = 80.1;
/// Distance of the contact-reaction probability [m] (Geant4's `Rs`).
pub const CONTACT_DISTANCE: f64 = 0.29e-9;
/// M^-1 s^-1 per m^3 mol^-1 s^-1.
const PER_MOLAR: f64 = 1e3;

/// Where the type of a pair reaction comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeSource {
    /// The `irt_type` key of the network.
    Network,
    /// The rule `k_obs >= k_dif`.
    Rule,
    /// A reactant is missing from the species table: fully
    /// diffusion-controlled.
    Fallback,
}

/// Reaction between two particles.
#[derive(Debug, Clone, PartialEq)]
pub struct PairReaction {
    /// Index of the reaction in the network (`SimMatrices::k`).
    pub index: usize,
    /// Reactant species, `a <= b`.
    pub a: usize,
    pub b: usize,
    pub kind: IrtType,
    pub source: TypeSource,
    /// Rate constant of the reaction events [M^-1 s^-1].
    pub k_obs: f64,
    /// Sum of the two diffusion coefficients [m^2/s].
    pub d_sum: f64,
    /// Reaction radius [m]: the sum of the van der Waals radii of a partially
    /// diffusion-controlled reaction, else the effective radius.
    pub r_reaction: f64,
    /// Onsager radius [m] (0 unless both reactants are charged; negative for
    /// opposite charges).
    pub r_onsager: f64,
    /// Effective reaction radius [m].
    pub r_eff: f64,
    /// Diffusion rate constant [M^-1 s^-1].
    pub k_dif: f64,
    /// Activation rate constant [M^-1 s^-1] (infinite when fully
    /// diffusion-controlled).
    pub k_act: f64,
    /// Probability that a pair found inside `r_eff` at the start reacts.
    pub p_contact: f64,
    /// Product species, one entry per molecule.
    pub products: Vec<usize>,
}

/// Reaction of one particle: first order, or with a fixed or background
/// species.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeReaction {
    /// Index of the reaction in the network.
    pub index: usize,
    pub species: usize,
    /// Rate [1/s].
    pub rate: f64,
    /// Background species consumed with the particle, if any.
    pub background: Option<usize>,
    /// Product species, one entry per molecule.
    pub products: Vec<usize>,
}

/// The network as seen by the IRT stage.
#[derive(Debug, Clone, PartialEq)]
pub struct IrtModel {
    /// Labels of the dynamic species, by index.
    pub species: Vec<String>,
    /// Diffusion coefficient of each species [m^2/s].
    pub diffusion: Vec<f64>,
    /// Background concentration of each species [mol/L].
    pub background: Vec<f64>,
    pub pairs: Vec<PairReaction>,
    pub molecule: Vec<MoleculeReaction>,
    /// Network reactions without a particle reactant, with the reason.
    pub skipped: Vec<(usize, String)>,
    /// Warnings raised while building (also logged).
    pub warnings: Vec<String>,
    /// Pair reactions of each species pair, `a * n + b` with `a <= b`.
    by_pair: Vec<Vec<u32>>,
    /// Molecule reactions of each species.
    by_species: Vec<Vec<u32>>,
}

/// Onsager radius [m] of two charges (Geant4: 293.15 K, permittivity 80.1).
pub fn onsager_radius(charge_a: i32, charge_b: i32) -> f64 {
    let q = f64::from(charge_a * charge_b);
    if q == 0.0 {
        return 0.0;
    }
    q * E_SI * E_SI
        / (4.0 * PI * EPSILON0 * K_BOLTZMANN)
        / (ONSAGER_TEMPERATURE * WATER_PERMITTIVITY)
}

/// Effective radius [m] of a reaction radius `r` and an Onsager radius `rc`.
pub fn effective_radius(r: f64, rc: f64) -> f64 {
    if rc == 0.0 {
        r
    } else {
        -rc / (1.0 - (rc / r).exp())
    }
}

/// Radii and rates of a pair reaction.
struct PairRates {
    kind: IrtType,
    source: TypeSource,
    r_reaction: f64,
    r_onsager: f64,
    r_eff: f64,
    k_dif: f64,
    k_act: f64,
    p_contact: f64,
}

impl PairRates {
    /// Fully diffusion-controlled: `Reff = k_obs / (4 pi D NA)`, with `D` the
    /// sum of the two coefficients, or the single one for identical reactants.
    fn full(k_obs: f64, d_sum: f64, identical: bool, source: TypeSource) -> Self {
        let d = if identical { 0.5 * d_sum } else { d_sum };
        let r_eff = k_obs / PER_MOLAR / (4.0 * PI * d * AVOGADRO);
        Self {
            kind: IrtType::Full,
            source,
            r_reaction: r_eff,
            r_onsager: 0.0,
            r_eff,
            k_dif: k_obs,
            k_act: f64::INFINITY,
            p_contact: 1.0,
        }
    }
}

/// `(r_reaction, r_onsager, r_eff, k_dif)` of a partially
/// diffusion-controlled pair.
fn partial_radii(
    radius_sum: f64,
    charges: (i32, i32),
    d_sum: f64,
    identical: bool,
) -> (f64, f64, f64, f64) {
    let rc = onsager_radius(charges.0, charges.1);
    let r_eff = effective_radius(radius_sum, rc);
    let mut k_dif = 4.0 * PI * d_sum * r_eff * AVOGADRO * PER_MOLAR;
    if identical {
        k_dif *= 0.5;
    }
    (radius_sum, rc, r_eff, k_dif)
}

impl IrtModel {
    /// Model of the network of `env` with the species properties of `table`.
    ///
    /// The diffusion coefficient of a species comes from `table`; a species
    /// absent from it does not diffuse and its pair reactions are fully
    /// diffusion-controlled (with a warning).
    pub fn new(env: &SimEnv, mat: &SimMatrices, table: &SpeciesTable) -> Result<Self> {
        if !env.reactions.enzymatic.is_empty() {
            return Err(RadkinError::Config(
                "irt: Michaelis-Menten reactions are not supported in the IRT stage".into(),
            ));
        }
        let n = mat.n_dyn;
        let species = env.species.sorted_dyn_labels();
        let props: Vec<_> = species.iter().map(|l| table.get(l).copied()).collect();
        let diffusion: Vec<f64> = props
            .iter()
            .map(|p| p.map_or(0.0, |p| p.diffusion))
            .collect();
        let mut background = vec![0.0; n];
        for s in env.species.dyn_sp.values() {
            background[s.index] = s.initial_cc(env.temperature);
        }
        let mut model = Self {
            species,
            diffusion,
            background,
            pairs: Vec::new(),
            molecule: Vec::new(),
            skipped: Vec::new(),
            warnings: Vec::new(),
            by_pair: vec![Vec::new(); n * n],
            by_species: vec![Vec::new(); n],
        };
        let mut missing: Vec<usize> = Vec::new();
        for (j, r) in env.reactions.iter_all_k().enumerate() {
            if !mat.k.negative_index[j].is_empty() {
                return Err(RadkinError::Config(format!(
                    "irt: reaction '{r}' has a negative product, which the IRT stage does not support"
                )));
            }
            let mut fixed = 1.0;
            let mut reactants: Vec<usize> = Vec::new();
            for &i in &mat.k.reactant_index[j] {
                if i < n {
                    reactants.push(i);
                } else if i < mat.n_species {
                    fixed *= mat.cst_sp[i - n];
                }
            }
            let mut products = Vec::new();
            for s in &r.products {
                if s.ope != Operation::Add || env.species.is_cst(&s.label) {
                    continue;
                }
                let i = env.species.index_of(&s.label)?;
                products.extend(std::iter::repeat_n(i, s.stoi as usize));
            }
            let k = mat.k.k_values[j] * fixed;
            match reactants.as_slice() {
                [] => model
                    .skipped
                    .push((j, format!("'{r}' has no particle reactant"))),
                [s] => model.push_molecule(MoleculeReaction {
                    index: j,
                    species: *s,
                    rate: k,
                    background: None,
                    products,
                }),
                [x, y] => {
                    let (a, b) = (*x.min(y), *x.max(y));
                    let identical = a == b;
                    let d_sum = model.diffusion[a] + model.diffusion[b];
                    if d_sum.is_nan() || d_sum <= 0.0 {
                        return Err(RadkinError::Config(format!(
                            "irt: reaction '{r}' has two reactants that do not diffuse (no diffusion coefficient)"
                        )));
                    }
                    let rates = match (props[a], props[b]) {
                        (Some(pa), Some(pb)) => {
                            let (r_reaction, r_onsager, r_eff, k_dif) = partial_radii(
                                pa.radius + pb.radius,
                                (pa.charge, pb.charge),
                                d_sum,
                                identical,
                            );
                            let limited = k >= k_dif;
                            let (kind, source) = match r.irt_type {
                                Some(t) => (t, TypeSource::Network),
                                None if limited => (IrtType::Full, TypeSource::Rule),
                                None => (IrtType::Partial, TypeSource::Rule),
                            };
                            if kind == IrtType::Partial && limited {
                                model.warn(format!(
                                    "irt: reaction '{r}' is typed partially diffusion-controlled but its rate constant {k:.4e} is not below its diffusion limit {k_dif:.4e} M-1 s-1: it runs fully diffusion-controlled"
                                ));
                                PairRates::full(k, d_sum, identical, source)
                            } else if kind == IrtType::Full {
                                PairRates::full(k, d_sum, identical, source)
                            } else {
                                let k_act = k_dif * k / (k_dif - k);
                                let rs = CONTACT_DISTANCE;
                                PairRates {
                                    kind,
                                    source,
                                    r_reaction,
                                    r_onsager,
                                    r_eff,
                                    k_dif,
                                    k_act,
                                    p_contact: rs / (rs + (k_dif / k_act) * (r_eff + rs)),
                                }
                            }
                        }
                        _ => {
                            for s in [a, b] {
                                if props[s].is_none() && !missing.contains(&s) {
                                    missing.push(s);
                                }
                            }
                            PairRates::full(k, d_sum, identical, TypeSource::Fallback)
                        }
                    };
                    let idx = model.pairs.len() as u32;
                    model.by_pair[a * n + b].push(idx);
                    model.pairs.push(PairReaction {
                        index: j,
                        a,
                        b,
                        kind: rates.kind,
                        source: rates.source,
                        k_obs: k,
                        d_sum,
                        r_reaction: rates.r_reaction,
                        r_onsager: rates.r_onsager,
                        r_eff: rates.r_eff,
                        k_dif: rates.k_dif,
                        k_act: rates.k_act,
                        p_contact: rates.p_contact,
                        products: products.clone(),
                    });
                    // A particle also meets the background of its partner.
                    let with_background = |species: usize, partner: usize, factor: f64| {
                        (model.background[partner] > 0.0).then(|| MoleculeReaction {
                            index: j,
                            species,
                            rate: factor * k * model.background[partner],
                            background: Some(partner),
                            products: products.clone(),
                        })
                    };
                    let extra: Vec<MoleculeReaction> = if identical {
                        with_background(a, a, 2.0).into_iter().collect()
                    } else {
                        with_background(a, b, 1.0)
                            .into_iter()
                            .chain(with_background(b, a, 1.0))
                            .collect()
                    };
                    for m in extra {
                        model.push_molecule(m);
                    }
                }
                _ => {
                    return Err(RadkinError::Config(format!(
                        "irt: reaction '{r}' has more than two particle reactants, which the IRT stage does not support"
                    )));
                }
            }
        }
        for s in missing {
            model.warn(format!(
                "irt: species '{}' is not in the species table (radius, charge): its reactions run fully diffusion-controlled",
                model.species[s]
            ));
        }
        Ok(model)
    }

    fn warn(&mut self, msg: String) {
        log::warn!("{msg}");
        self.warnings.push(msg);
    }

    fn push_molecule(&mut self, m: MoleculeReaction) {
        self.by_species[m.species].push(self.molecule.len() as u32);
        self.molecule.push(m);
    }

    pub fn n_species(&self) -> usize {
        self.species.len()
    }

    /// Pair reactions (indices in `pairs`) of the species `a` and `b`.
    pub fn pairs_of(&self, a: usize, b: usize) -> &[u32] {
        let (a, b) = (a.min(b), a.max(b));
        &self.by_pair[a * self.n_species() + b]
    }

    /// Molecule reactions (indices in `molecule`) of a species.
    pub fn molecule_of(&self, species: usize) -> &[u32] {
        &self.by_species[species]
    }

    /// Largest effective radius of the pair reactions [m].
    pub fn max_radius(&self) -> f64 {
        self.pairs.iter().map(|p| p.r_eff).fold(0.0, f64::max)
    }

    /// Largest sum of two diffusion coefficients of a reacting pair [m^2/s].
    pub fn max_d_sum(&self) -> f64 {
        self.pairs.iter().map(|p| p.d_sum).fold(0.0, f64::max)
    }
}
