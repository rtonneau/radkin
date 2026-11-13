//! Reaction objects (port of radiopyo `basic_reaction.py`, `k_reactions.py`,
//! `acid_base.py`, `e_reactions.py`, `michaelis.py`, `radiolytic.py`).

use std::fmt;

use crate::chemistry::parser::{GenericReaction, Operation, Side, revert};
use crate::consts::{AVOGADRO, ELEMENTARY_CHARGE, GAS_CONSTANT, ZERO_CELSIUS};
use crate::error::{RadkinError, Result};
use crate::poly::Polynomial;

/// A species of a reaction side.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionSpecies {
    pub label: String,
    pub stoi: u32,
    pub ope: Operation,
}

fn side(s: &Side) -> Vec<ReactionSpecies> {
    s.iter()
        .map(|(label, (stoi, ope))| ReactionSpecies {
            label: label.clone(),
            stoi: *stoi,
            ope: *ope,
        })
        .collect()
}

/// radiopyo `ReactionSpeciesCollection.as_string`.
pub fn side_string(species: &[ReactionSpecies]) -> String {
    let mut out: Vec<String> = Vec::new();
    for sp in species {
        let stoi = if sp.stoi <= 1 {
            String::new()
        } else {
            sp.stoi.to_string()
        };
        if !out.is_empty() || sp.ope == Operation::Subtract {
            out.push(sp.ope.to_string());
        }
        out.push(format!("{stoi}{}", sp.label));
    }
    out.join(" ")
}

/// `ge` [#/100 eV] to a radiolytic rate constant [mol/L/Gy].
pub fn ge_to_kr(ge: f64) -> f64 {
    ge * 1.0 / ELEMENTARY_CHARGE / 100.0 / AVOGADRO
}

/// Inverse of [`ge_to_kr`].
pub fn kr_to_ge(kr: f64) -> f64 {
    kr / 1.0 * ELEMENTARY_CHARGE * 100.0 * AVOGADRO
}

/// Rate constant of a k-reaction.
#[derive(Debug, Clone, PartialEq)]
pub enum KValue {
    /// `log10(k)` as a polynomial of `1/T` [1/K].
    Log10Poly(Polynomial),
    /// `k = A exp(-Ea * 1000 / (R T))`, `Ea` in kJ/mol.
    Arrhenius { a: f64, ea: f64 },
}

impl KValue {
    pub fn at(&self, t_celsius: f64) -> f64 {
        let kelvin = t_celsius + ZERO_CELSIUS;
        match self {
            Self::Log10Poly(p) => 10f64.powf(p.eval(1.0 / kelvin)),
            Self::Arrhenius { a, ea } => a * (-ea * 1000.0 / (GAS_CONSTANT * kelvin)).exp(),
        }
    }
}

/// Flavour of a k-reaction (label prefix and display).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KKind {
    K,
    Arrhenius,
    Constant,
    /// Forward or backward reaction of an acid/base or equilibrium pair.
    Pair,
}

/// Mass-action reaction.
#[derive(Debug, Clone, PartialEq)]
pub struct KReaction {
    pub reactants: Vec<ReactionSpecies>,
    pub products: Vec<ReactionSpecies>,
    pub k: Option<KValue>,
    pub kind: KKind,
    /// Global reaction index.
    pub index: usize,
}

impl KReaction {
    fn from_generic(g: &GenericReaction, k: Option<KValue>, kind: KKind) -> Self {
        let reactants = if kind == KKind::Constant {
            Vec::new()
        } else {
            side(&g.reactants)
        };
        Self {
            reactants,
            products: side(&g.products),
            k,
            kind,
            index: usize::MAX,
        }
    }

    /// k-reaction (`[[reactions.k_reaction]]`) with the n·k rule (ADR 0001).
    pub fn new_k(g: &GenericReaction, k_value: Option<Polynomial>) -> Self {
        let mut r = Self::from_generic(g, k_value.map(KValue::Log10Poly), KKind::K);
        r.apply_nk_convention();
        r
    }

    /// Arrhenius reaction.
    pub fn new_arrhenius(g: &GenericReaction, a: f64, ea: f64) -> Self {
        let r = Self::from_generic(g, Some(KValue::Arrhenius { a, ea }), KKind::Arrhenius);
        r.warn_nk_not_divided();
        r
    }

    /// Constant production `-> X` (no reactant).
    pub fn new_constant(g: &GenericReaction, k_value: Option<Polynomial>) -> Self {
        Self::from_generic(g, k_value.map(KValue::Log10Poly), KKind::Constant)
    }

    fn new_pair(g: &GenericReaction, k: Polynomial) -> Self {
        Self::from_generic(g, Some(KValue::Log10Poly(k)), KKind::Pair)
    }

    fn warn_nk_not_divided(&self) {
        if self.reactants.len() == 1 && self.reactants[0].stoi > 1 {
            log::warn!(
                "Reaction '{self}': single reactant with stoichiometry {} and a non-scalar k (polynomial, Arrhenius or equilibrium): k is used as given, not divided by n (ADR 0001).",
                self.reactants[0].stoi
            );
        }
    }

    fn apply_nk_convention(&mut self) {
        let scalar = matches!(&self.k, Some(KValue::Log10Poly(p)) if p.len() == 1);
        if scalar && self.reactants.len() == 1 {
            let n = self.reactants[0].stoi as f64;
            let k = self
                .k
                .as_ref()
                .unwrap()
                .at(crate::consts::DEFAULT_TEMPERATURE)
                / n;
            self.k = Some(KValue::Log10Poly(Polynomial::from_value(k, true)));
        } else if self.k.is_some() {
            self.warn_nk_not_divided();
        }
    }

    /// Rate constant at `t_celsius`.
    pub fn k_value(&self, t_celsius: f64) -> Result<f64> {
        self.k.as_ref().map(|k| k.at(t_celsius)).ok_or_else(|| {
            RadkinError::Config(format!(
                "Unable to compute reaction constant. No data provided for reaction: {self}"
            ))
        })
    }

    /// Sum of the reactant stoichiometries.
    pub fn tot_number_reactants(&self) -> u32 {
        self.reactants.iter().map(|s| s.stoi).sum()
    }

    pub fn as_label(&self) -> String {
        let prefix = if self.kind == KKind::Constant {
            "CR"
        } else {
            "KR"
        };
        format!("{prefix}_{}", self.index)
    }

    pub fn has_reactant(&self, label: &str) -> bool {
        self.reactants.iter().any(|s| s.label == label)
    }

    pub fn has_product(&self, label: &str) -> bool {
        self.products.iter().any(|s| s.label == label)
    }

    pub fn contains(&self, label: &str) -> bool {
        self.has_reactant(label) || self.has_product(label)
    }

    /// Stoichiometry of `label` (reactant first, then product).
    pub fn stoi_of(&self, label: &str) -> Option<u32> {
        self.reactants
            .iter()
            .chain(&self.products)
            .find(|s| s.label == label)
            .map(|s| s.stoi)
    }
}

impl fmt::Display for KReaction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let products = side_string(&self.products);
        match self.kind {
            KKind::Constant => write!(f, "-> {products}"),
            KKind::Pair => write!(f, "{} <-> {products}", side_string(&self.reactants)),
            _ => write!(f, "{} -> {products}", side_string(&self.reactants)),
        }
    }
}

/// Acid/base data of a reaction pair.
#[derive(Debug, Clone, PartialEq)]
pub struct AcidBaseCouple {
    /// pKa polynomial of the temperature [degC].
    pub pka: Polynomial,
    pub acid: String,
    pub base: String,
}

/// Forward/backward pair: acid/base (`acid_base`) or equilibrium.
#[derive(Debug, Clone, PartialEq)]
pub struct BiReaction {
    pub forward: KReaction,
    pub backward: KReaction,
    pub acid_base: Option<AcidBaseCouple>,
    /// `log10(K_eq)` polynomial of `1/T` (equilibrium only, when defined).
    pub k_equilibrium: Option<Polynomial>,
}

impl BiReaction {
    fn pair(g: &GenericReaction, kf: Polynomial, kb: Polynomial) -> Result<(KReaction, KReaction)> {
        let forward = KReaction::new_pair(g, kf);
        let backward = KReaction::new_pair(&revert(g)?, kb);
        forward.warn_nk_not_divided();
        backward.warn_nk_not_divided();
        Ok((forward, backward))
    }

    /// radiopyo `AcidBase`.
    pub fn acid_base(
        g: &GenericReaction,
        pka: Polynomial,
        kf: Polynomial,
        kb: Polynomial,
    ) -> Result<Self> {
        let (forward, backward) = Self::pair(g, kf, kb)?;
        let (acid, base) = couple(&forward, &backward)?;
        Ok(Self {
            forward,
            backward,
            acid_base: Some(AcidBaseCouple { pka, acid, base }),
            k_equilibrium: None,
        })
    }

    /// radiopyo `EReaction`: two of k_forward, k_backward, k_equilibrium.
    pub fn equilibrium(
        g: &GenericReaction,
        kf: Option<Polynomial>,
        kb: Option<Polynomial>,
        keq: Option<Polynomial>,
    ) -> Result<Self> {
        let missing = || {
            RadkinError::Config(format!(
                "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward): {}",
                side_string(&side(&g.reactants))
            ))
        };
        let (kf, kb, keq) = match keq {
            None => {
                let (kf, kb) = (kf.ok_or_else(missing)?, kb.ok_or_else(missing)?);
                let keq = (kf.coef[0].is_finite() && kb.coef[0].is_finite()).then(|| kf.sub(&kb));
                (kf, kb, keq)
            }
            Some(keq) => match (kf, kb) {
                (_, Some(kb)) => (keq.add(&kb), kb, Some(keq)),
                (Some(kf), None) => {
                    let kb = kf.sub(&keq);
                    (kf, kb, Some(keq))
                }
                (None, None) => return Err(missing()),
            },
        };
        let (forward, backward) = Self::pair(g, kf, kb)?;
        Ok(Self {
            forward,
            backward,
            acid_base: None,
            k_equilibrium: keq,
        })
    }
}

fn couple(forward: &KReaction, backward: &KReaction) -> Result<(String, String)> {
    let is_acid = forward.contains("H3O_plus") || forward.contains("H_plus");
    let is_base = forward.contains("OH_minus");
    let member = |r: &KReaction,
                  coll: &[ReactionSpecies],
                  kind: &str,
                  partners: &[&str]|
     -> Result<String> {
        if coll.len() > 2 {
            return Err(RadkinError::Config(format!(
                "Unable to identify acid/base couple of {r}. The reaction has more than 2 different {kind}s"
            )));
        }
        if coll.len() == 2 && !coll.iter().any(|s| partners.contains(&s.label.as_str())) {
            return Err(RadkinError::Config(format!(
                "Unable to identify acid/base couple of {r}. The reaction has more than 1 {kind} but none of them is {}.",
                partners.join(" or ")
            )));
        }
        let members: Vec<&ReactionSpecies> = coll
            .iter()
            .filter(|s| !partners.contains(&s.label.as_str()))
            .collect();
        Ok(members
            .last()
            .copied()
            .or(coll.last())
            .map(|s| s.label.clone())
            .unwrap_or_default())
    };
    if is_acid {
        let r = if forward.has_product("H3O_plus") || forward.has_product("H_plus") {
            forward
        } else {
            backward
        };
        let acid = member(r, &r.reactants, "reactant", &["H2O"])?;
        let base = member(r, &r.products, "product", &["H3O_plus", "H_plus"])?;
        return Ok((acid, base));
    }
    if is_base {
        let r = if forward.has_reactant("OH_minus") {
            forward
        } else {
            backward
        };
        let base = member(r, &r.products, "product", &["H2O"])?;
        let acid = member(r, &r.reactants, "reactant", &["OH_minus"])?;
        return Ok((acid, base));
    }
    Err(RadkinError::Config(format!(
        "Unable to identify Acid or Base nature of {forward}"
    )))
}

/// Michaelis–Menten enzymatic reaction.
#[derive(Debug, Clone, PartialEq)]
pub struct MichaelisMenten {
    pub reactants: Vec<ReactionSpecies>,
    pub products: Vec<ReactionSpecies>,
    pub enzyme: String,
    /// k_cat [1/s] (used as given, ADR 0006).
    pub k_value: f64,
    /// Michaelis constant [mol/L].
    pub k_micha: f64,
    pub index: usize,
}

impl MichaelisMenten {
    pub fn new(g: &GenericReaction, k_value: f64, k_micha: f64) -> Result<Self> {
        let enzyme = g
            .enzyme
            .clone()
            .ok_or_else(|| RadkinError::Config("No enzyme found for enzymatic reaction".into()))?;
        if g.reactants.len() > 1 {
            return Err(RadkinError::Config(
                "MichaelisMenten reaction must have only one substrate".into(),
            ));
        }
        Ok(Self {
            reactants: side(&g.reactants),
            products: side(&g.products),
            enzyme,
            k_value,
            k_micha,
            index: usize::MAX,
        })
    }

    pub fn substrate(&self) -> &str {
        &self.reactants[0].label
    }

    pub fn as_label(&self) -> String {
        format!("EZ_{}", self.index)
    }

    pub fn has_reactant(&self, label: &str) -> bool {
        self.reactants.iter().any(|s| s.label == label)
    }

    pub fn has_product(&self, label: &str) -> bool {
        self.products.iter().any(|s| s.label == label)
    }

    pub fn stoi_of(&self, label: &str) -> Option<u32> {
        self.reactants
            .iter()
            .chain(&self.products)
            .find(|s| s.label == label)
            .map(|s| s.stoi)
    }
}

impl fmt::Display for MichaelisMenten {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} -- {} >> {}",
            side_string(&self.reactants),
            self.enzyme,
            side_string(&self.products)
        )
    }
}

/// Production of a species by the radiation.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiolyticReaction {
    pub species: String,
    /// Rate constant [mol/L/Gy].
    pub kr: f64,
    pub index: usize,
}

impl RadiolyticReaction {
    pub fn from_ge(species: &str, ge: f64) -> Self {
        Self {
            species: species.to_string(),
            kr: ge_to_kr(ge),
            index: usize::MAX,
        }
    }

    pub fn ge(&self) -> f64 {
        kr_to_ge(self.kr)
    }

    pub fn update_from_ge(&mut self, ge: f64) {
        self.kr = ge_to_kr(ge);
    }

    /// `G_<species>`.
    pub fn label(&self) -> String {
        format!("G_{}", self.species)
    }

    /// Concentration produced over `t0` at `dose_rate` [mol/L].
    pub fn initial_cc_before(&self, t0: f64, dose_rate: f64) -> f64 {
        self.kr * dose_rate * t0
    }
}

impl fmt::Display for RadiolyticReaction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RadiolyticReaction('{}')", self.species)
    }
}
