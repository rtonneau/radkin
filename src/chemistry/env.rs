//! Simulation environment: species, reactions and parameters (port of
//! radiopyo `simulation/environment_build.py`).

use crate::chemistry::parser::{parse_reaction, parse_species};
use crate::chemistry::reactions::{BiReaction, KReaction, MichaelisMenten, RadiolyticReaction};
use crate::chemistry::species::ListSpecies;
use crate::config::{Config, Section};
use crate::consts::H2O_CC;
use crate::error::{RadkinError, Result};
use crate::gvalues::GValuesBoscolo;
use crate::poly::Polynomial;

/// Every reaction, with global indices: acid/base (2 each), equilibrium
/// (2 each), k-reactions, enzymatic, radiolytic.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReactionsList {
    pub acid_base: Vec<BiReaction>,
    pub equilibrium: Vec<BiReaction>,
    /// `k_reaction`, then `arrhenius`, then `constant` entries.
    pub k: Vec<KReaction>,
    pub enzymatic: Vec<MichaelisMenten>,
    pub radiolytic: Vec<RadiolyticReaction>,
}

impl ReactionsList {
    fn assign_indices(&mut self) {
        let mut idx = 0;
        for bi in self.acid_base.iter_mut().chain(self.equilibrium.iter_mut()) {
            bi.forward.index = idx;
            bi.backward.index = idx + 1;
            idx += 2;
        }
        for r in &mut self.k {
            r.index = idx;
            idx += 1;
        }
        for r in &mut self.enzymatic {
            r.index = idx;
            idx += 1;
        }
        for r in &mut self.radiolytic {
            r.index = idx;
            idx += 1;
        }
    }

    pub fn total_kreactions(&self) -> usize {
        self.k.len() + 2 * (self.acid_base.len() + self.equilibrium.len())
    }

    /// Pairs (acid/base, then equilibrium).
    pub fn iter_pairs(&self) -> impl Iterator<Item = &BiReaction> {
        self.acid_base.iter().chain(&self.equilibrium)
    }

    /// Every mass-action reaction in index order (pairs fwd/bwd, then k).
    pub fn iter_all_k(&self) -> impl Iterator<Item = &KReaction> {
        self.iter_pairs()
            .flat_map(|b| [&b.forward, &b.backward])
            .chain(&self.k)
    }
}

/// The chemistry of a simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct SimEnv {
    pub species: ListSpecies,
    pub reactions: ReactionsList,
    /// Temperature [degC], the single source of truth.
    pub temperature: f64,
    pub ph: f64,
    pub force_cst_ph: bool,
    pub variable_gvalues: bool,
}

/// Species of the pH.
pub const PH_SPECIES: &str = "H3O_plus";

fn number(v: &toml::Value, what: &str) -> Result<f64> {
    match v {
        toml::Value::Float(f) => Ok(*f),
        toml::Value::Integer(i) => Ok(*i as f64),
        _ => Err(RadkinError::Config(format!(
            "{what}: expected a number, got {v}"
        ))),
    }
}

fn numbers(v: &toml::Value, what: &str) -> Result<Vec<f64>> {
    match v {
        toml::Value::Array(a) => a.iter().map(|x| number(x, what)).collect(),
        _ => Ok(vec![number(v, what)?]),
    }
}

/// radiopyo `polynomial_factory` from a config value.
fn poly(v: &toml::Value, log10: bool, what: &str) -> Result<Polynomial> {
    match v {
        toml::Value::Array(_) => {
            let c = numbers(v, what)?;
            if c.is_empty() {
                return Err(RadkinError::Config(format!(
                    "{what}: cannot make a polynomial with no coeff"
                )));
            }
            Ok(Polynomial::from_coefs(c, log10))
        }
        _ => Ok(Polynomial::from_value(number(v, what)?, log10)),
    }
}

fn reaction_str(r: &toml::Table) -> Result<&str> {
    r.get("reaction")
        .and_then(|v| v.as_str())
        .ok_or_else(|| RadkinError::Config("reaction entry without a 'reaction' string".into()))
}

fn opt_poly(r: &toml::Table, key: &str, log10: bool) -> Result<Option<Polynomial>> {
    r.get(key).map(|v| poly(v, log10, key)).transpose()
}

fn require<'a>(r: &'a toml::Table, key: &str, what: &str) -> Result<&'a toml::Value> {
    r.get(key)
        .ok_or_else(|| RadkinError::Config(format!("{what}: missing '{key}' ({})", r)))
}

fn add_species_of(sp: &mut ListSpecies, labels: impl IntoIterator<Item = String>) {
    for l in labels {
        sp.get_or_create_dyn(&l);
    }
}

impl SimEnv {
    /// Build the environment from a resolved configuration (radiopyo
    /// `SimEnv._chemistry_builder`).
    pub fn from_config(cfg: &Config) -> Result<Self> {
        let ph = cfg.ph();
        let force_cst_ph = cfg.force_cst_ph();
        let mut sp = ListSpecies::new();
        for (label, cc) in &cfg.fixed {
            sp.add_cst(&parse_species(label)?, poly(cc, false, label)?)?;
        }
        sp.add_cst("H2O", Polynomial::from_coefs(H2O_CC.to_vec(), false))?;
        let h3o = Polynomial::from_value(10f64.powf(-ph), false);
        let oh = Polynomial::from_value(10f64.powf(-14.0 + ph), false);
        if force_cst_ph {
            if !sp.contains(PH_SPECIES) {
                sp.add_cst(PH_SPECIES, h3o)?;
            }
            if !sp.contains("OH_minus") {
                sp.add_cst("OH_minus", oh)?;
            }
        } else {
            sp.add_dyn(PH_SPECIES, Some(h3o))?;
            sp.add_dyn("OH_minus", Some(oh))?;
        }

        let mut reactions = ReactionsList::default();
        let rx = |s: Section| cfg.reactions[&s].iter();
        for r in rx(Section::AcidBase) {
            let g = parse_reaction(reaction_str(r)?)?;
            if r.contains_key("k_value") {
                return Err(RadkinError::Config(
                    "'k_value' is ambiguous for an AcidBase reaction. Use k_forward or k_backward instead.".into(),
                ));
            }
            let pka = Polynomial::from_coefs(
                numbers(
                    require(
                        r,
                        "pKa",
                        "Need a reaction string AND a pKa to create an Acid Base rel.",
                    )?,
                    "pKa",
                )?,
                false,
            );
            let kf = opt_poly(r, "k_forward", true)?.ok_or_else(|| {
                RadkinError::Config("For A/B reaction, k_forward must be defined".into())
            })?;
            let kb = opt_poly(r, "k_backward", true)?.ok_or_else(|| {
                RadkinError::Config("For A/B reaction, k_backward must be defined".into())
            })?;
            add_species_of(
                &mut sp,
                g.reactants.keys().chain(g.products.keys()).cloned(),
            );
            reactions
                .acid_base
                .push(BiReaction::acid_base(&g, pka, kf, kb)?);
        }
        for r in rx(Section::Equilibrium) {
            if r.contains_key("k_value") {
                return Err(RadkinError::Config(
                    "'k_value' is ambiguous for an Equilibrium Reaction. Use k_forward or k_backward instead.".into(),
                ));
            }
            let g = parse_reaction(reaction_str(r)?)?;
            add_species_of(
                &mut sp,
                g.reactants.keys().chain(g.products.keys()).cloned(),
            );
            reactions.equilibrium.push(BiReaction::equilibrium(
                &g,
                opt_poly(r, "k_forward", true)?,
                opt_poly(r, "k_backward", true)?,
                opt_poly(r, "k_equilibrium", true)?,
            )?);
        }
        // Species are created in the order k_reaction, constant, arrhenius;
        // reactions are listed k_reaction, arrhenius, constant (radiopyo).
        let mut k_list = Vec::new();
        for r in rx(Section::KReaction) {
            let g = parse_reaction(reaction_str(r)?)?;
            let k = poly(require(r, "k_value", reaction_str(r)?)?, true, "k_value")?;
            add_species_of(
                &mut sp,
                g.reactants.keys().chain(g.products.keys()).cloned(),
            );
            k_list.push(KReaction::new_k(&g, Some(k)));
        }
        let mut cst_list = Vec::new();
        for r in rx(Section::Constant) {
            let g = parse_reaction(reaction_str(r)?)?;
            let k = poly(require(r, "k_value", reaction_str(r)?)?, true, "k_value")?;
            add_species_of(&mut sp, g.products.keys().cloned());
            cst_list.push(KReaction::new_constant(&g, Some(k)));
        }
        let mut arr_list = Vec::new();
        for r in rx(Section::Arrhenius) {
            let g = parse_reaction(reaction_str(r)?)?;
            let a = number(require(r, "A", reaction_str(r)?)?, "A")?;
            let ea = number(require(r, "Ea", reaction_str(r)?)?, "Ea")?;
            add_species_of(
                &mut sp,
                g.reactants.keys().chain(g.products.keys()).cloned(),
            );
            arr_list.push(KReaction::new_arrhenius(&g, a, ea));
        }
        reactions.k = k_list.into_iter().chain(arr_list).chain(cst_list).collect();
        for r in rx(Section::Enzymatic) {
            let g = parse_reaction(reaction_str(r)?)?;
            let k_value = number(require(r, "k_value", reaction_str(r)?)?, "k_value")?;
            let k_micha = number(require(r, "k_micha", reaction_str(r)?)?, "k_micha")?;
            let ez = MichaelisMenten::new(&g, k_value, k_micha)?;
            sp.get_or_create_dyn(&ez.enzyme);
            add_species_of(
                &mut sp,
                g.reactants.keys().chain(g.products.keys()).cloned(),
            );
            reactions.enzymatic.push(ez);
        }
        for (product, ge) in &cfg.radiolytic {
            let label = parse_species(product)?;
            sp.get_or_create_dyn(&label);
            reactions
                .radiolytic
                .push(RadiolyticReaction::from_ge(&label, number(ge, product)?));
        }
        for (label, cc) in &cfg.initial {
            if let Some(s) = sp.get_mut(label.trim()) {
                s.cc_init = Some(poly(cc, false, label)?);
            }
        }
        if sp.is_dyn("trash") {
            sp.remove("trash");
            sp.add_cst("trash", Polynomial::from_value(0.0, false))?;
        }
        match cfg.species_order() {
            Some((d, c)) => sp.set_order(&d, &c)?,
            None => {
                let starter: Vec<String> = reactions
                    .acid_base
                    .iter()
                    .flat_map(|r| {
                        let c = r.acid_base.as_ref().expect("acid/base pair");
                        [c.acid.clone(), c.base.clone()]
                    })
                    .collect();
                sp.sort(&starter)?;
            }
        }
        reactions.assign_indices();
        Ok(Self {
            species: sp,
            reactions,
            temperature: cfg.temperature(),
            ph,
            force_cst_ph,
            variable_gvalues: cfg.variable_gvalues(),
        })
    }

    /// Initial concentrations of the dynamic species [µmol/L], by index.
    pub fn initial_values(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.species.n_dyn()];
        for s in self.species.dyn_sp.values() {
            out[s.index] = s.initial_cc(self.temperature) * 1e6;
        }
        out
    }

    /// Set the initial concentration [mol/L] of a species.
    pub fn set_initial_cc(&mut self, label: &str, cc: f64) -> Result<()> {
        let s = self
            .species
            .get_mut(label)
            .ok_or_else(|| RadkinError::Value(format!("unknown species '{label}'")))?;
        s.cc_init = Some(Polynomial::from_value(cc, false));
        Ok(())
    }

    /// Rate constants of every mass-action reaction, in index order.
    pub fn k_values(&self) -> Result<Vec<f64>> {
        self.reactions
            .iter_all_k()
            .map(|r| r.k_value(self.temperature))
            .collect()
    }

    /// `str(reaction)` of every reaction: k (pairs first), radiolytic, enzymatic.
    pub fn reaction_strings(&self) -> Vec<String> {
        let r = &self.reactions;
        r.iter_all_k()
            .map(|x| x.to_string())
            .chain(r.radiolytic.iter().map(|x| x.to_string()))
            .chain(r.enzymatic.iter().map(|x| x.to_string()))
            .collect()
    }

    /// `as_label()` of every reaction, same order as [`SimEnv::reaction_strings`].
    pub fn reaction_labels(&self) -> Vec<String> {
        let r = &self.reactions;
        r.iter_all_k()
            .map(|x| x.as_label())
            .chain(r.radiolytic.iter().map(|x| x.label()))
            .chain(r.enzymatic.iter().map(|x| x.as_label()))
            .collect()
    }

    /// Rewrite the G-values of the Boscolo species for `let_` [keV/µm]
    /// (radiopyo `UnitCell.change_LET`); `dose_rate_peak` scales them when
    /// `variable_GValues` is set.
    pub fn apply_let(
        &mut self,
        gv: &GValuesBoscolo,
        let_: f64,
        dose_rate_peak: Option<f64>,
    ) -> Result<()> {
        let dr = if self.variable_gvalues {
            dose_rate_peak
        } else {
            None
        };
        for r in &mut self.reactions.radiolytic {
            if !gv.contains(&r.species) || !self.species.is_dyn(&r.species) {
                continue;
            }
            r.update_from_ge(gv.value_at(&r.species, let_, dr)?);
        }
        Ok(())
    }
}
