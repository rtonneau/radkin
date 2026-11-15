//! Propensity model: the reaction network of [`SimMatrices`] on integer
//! molecule counts in a volume `V`.
//!
//! The LSODA path integrates rates `r` [mol/L/s] of concentrations; here a
//! reaction fires `r * N_A * V` times per second on average. For a
//! mass-action reaction with rate constant `k`, fixed reactants of
//! concentration `c_f` and `m` dynamic reactant slots, the propensity is
//!
//! ```text
//! a = k * prod(c_f) * N_A V * prod(falling factorial of the counts) / (N_A V)^m
//! ```
//!
//! so `2 A -> ...` (rate `k [A]^2`, stoichiometry -2, radiopyo's convention)
//! has `a = k / (N_A V) * n (n - 1)` and removes two molecules per firing.

use crate::consts::AVOGADRO;
use crate::error::{RadkinError, Result};
use crate::matrices::{ReactionPair, SimMatrices};

/// A reactant of a Michaelis–Menten reaction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Slot {
    /// Dynamic species index.
    Dyn(usize),
    /// Fixed species concentration [mol/L].
    Fixed(f64),
}

/// Rate law of a reaction.
#[derive(Debug, Clone, PartialEq)]
pub enum Law {
    /// `c * prod(falling factorial(n_s, mult))`.
    MassAction {
        /// `k * prod(fixed concentrations) * (N_A V)^(1 - m)` [1/s].
        c: f64,
        /// Dynamic reactants `(species, multiplicity)`.
        reactants: Vec<(usize, u32)>,
    },
    /// `k * E * S / (Km + S) * N_A V`, concentrations from the counts.
    Michaelis {
        k: f64,
        km: f64,
        enzyme: Slot,
        substrate: Slot,
    },
}

/// One reaction channel.
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    pub law: Law,
    /// Net change of every dynamic species it changes.
    pub delta: Vec<(usize, i64)>,
    /// `(species, count)` it needs to fire: every species it consumes on
    /// balance, "negative" products included.
    pub need: Vec<(usize, i64)>,
    /// Number of dynamic reactant slots (order seen by tau-leaping).
    pub order: u32,
}

/// Radiolytic production of one dynamic species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source {
    pub species: usize,
    /// Molecules per gray: `yield * N_A V`.
    pub per_gray: f64,
}

/// The network on molecule counts in one volume.
#[derive(Debug, Clone, PartialEq)]
pub struct SsaModel {
    pub n_dyn: usize,
    /// Volume [L].
    pub volume: f64,
    /// `N_A * V`: molecules per mol/L.
    pub nav: f64,
    /// Mass-action reactions (in `mat.k` order), then the enzymatic ones.
    pub reactions: Vec<Reaction>,
    pub sources: Vec<Source>,
    /// Reactions whose propensity reads the count of each species.
    pub readers: Vec<Vec<usize>>,
    /// Reactions to re-evaluate after reaction `j` fires.
    pub affects: Vec<Vec<usize>>,
    /// Reactions to re-evaluate after source `i` fires.
    pub source_affects: Vec<Vec<usize>>,
    /// Reversible pairs, declared with `<->` or inferred from one-way
    /// reactions ([`SimMatrices::reversible_pairs`]).
    pub pairs: Vec<ReactionPair>,
}

fn integer(coef: f64, what: &str) -> Result<i64> {
    if coef.fract() != 0.0 || coef.abs() > 1e9 {
        return Err(RadkinError::Value(format!(
            "ssa: {what} has a non-integer stoichiometry ({coef})"
        )));
    }
    Ok(coef as i64)
}

/// `(species, amount)` entries.
type Entries = Vec<(usize, i64)>;

/// Net change and needed counts of a stoichiometry column.
fn net(col: &[(usize, f64)], what: &str) -> Result<(Entries, Entries)> {
    let mut delta = Vec::with_capacity(col.len());
    let mut need = Vec::new();
    for &(s, c) in col {
        let d = integer(c, what)?;
        if d == 0 {
            continue;
        }
        delta.push((s, d));
        if d < 0 {
            need.push((s, -d));
        }
    }
    Ok((delta, need))
}

impl SsaModel {
    /// Model of `mat` in a volume of `volume` litres.
    pub fn new(mat: &SimMatrices, volume: f64) -> Result<Self> {
        if !(volume.is_finite() && volume > 0.0) {
            return Err(RadkinError::Value(format!(
                "ssa: volume must be finite and > 0 L, got {volume}"
            )));
        }
        let nav = AVOGADRO * volume;
        let n_dyn = mat.n_dyn;
        let slot = |i: usize| {
            if i < n_dyn {
                Slot::Dyn(i)
            } else {
                Slot::Fixed(mat.cst_sp[i - n_dyn])
            }
        };
        let mut reactions = Vec::new();
        for r in 0..mat.k.k_values.len() {
            let what = format!("mass-action reaction {r}");
            let mut c = mat.k.k_values[r] * nav;
            let mut reactants: Vec<(usize, u32)> = Vec::new();
            let mut order = 0;
            for &i in &mat.k.reactant_index[r] {
                match slot(i) {
                    Slot::Fixed(cc) => c *= cc,
                    Slot::Dyn(s) => {
                        c /= nav;
                        order += 1;
                        match reactants.iter_mut().find(|(sp, _)| *sp == s) {
                            Some(e) => e.1 += 1,
                            None => reactants.push((s, 1)),
                        }
                    }
                }
            }
            let (delta, need) = net(&mat.k.stoichiometry[r], &what)?;
            reactions.push(Reaction {
                law: Law::MassAction { c, reactants },
                delta,
                need,
                order,
            });
        }
        let ez = &mat.enzymatic;
        for r in 0..ez.k_values.len() {
            let what = format!("enzymatic reaction {r}");
            let (enzyme, substrate) = (slot(ez.enzyme_index[r]), slot(ez.substrate_index[r]));
            let order = [enzyme, substrate]
                .iter()
                .filter(|s| matches!(s, Slot::Dyn(_)))
                .count() as u32;
            let (delta, need) = net(&ez.stoichiometry[r], &what)?;
            reactions.push(Reaction {
                law: Law::Michaelis {
                    k: ez.k_values[r],
                    km: ez.k_michaelis[r],
                    enzyme,
                    substrate,
                },
                delta,
                need,
                order,
            });
        }
        let sources: Vec<Source> = mat
            .radiolytic
            .yields
            .iter()
            .zip(&mat.radiolytic.species)
            .filter_map(|(y, sp)| {
                sp.map(|species| Source {
                    species,
                    per_gray: y * nav,
                })
            })
            .collect();

        let mut readers = vec![Vec::new(); n_dyn];
        for (j, r) in reactions.iter().enumerate() {
            let mut read: Vec<usize> = r.need.iter().map(|(s, _)| *s).collect();
            match &r.law {
                Law::MassAction { reactants, .. } => read.extend(reactants.iter().map(|(s, _)| *s)),
                Law::Michaelis {
                    enzyme, substrate, ..
                } => {
                    for s in [enzyme, substrate] {
                        if let Slot::Dyn(i) = s {
                            read.push(*i);
                        }
                    }
                }
            }
            read.sort_unstable();
            read.dedup();
            for s in read {
                readers[s].push(j);
            }
        }
        let union = |species: &mut dyn Iterator<Item = usize>| {
            let mut v: Vec<usize> = species.flat_map(|s| readers[s].iter().copied()).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        let affects = reactions
            .iter()
            .map(|r| union(&mut r.delta.iter().map(|(s, _)| *s)))
            .collect();
        let source_affects = sources
            .iter()
            .map(|s| union(&mut std::iter::once(s.species)))
            .collect();
        Ok(Self {
            n_dyn,
            volume,
            nav,
            reactions,
            sources,
            readers,
            affects,
            source_affects,
            pairs: mat.reversible_pairs.clone(),
        })
    }

    /// Whether reaction `j` is bimolecular once the fixed species are
    /// dropped: mass action with two dynamic reactant slots.
    pub fn is_bimolecular(&self, j: usize) -> bool {
        let r = &self.reactions[j];
        matches!(r.law, Law::MassAction { .. }) && r.order == 2
    }

    /// Propensity [1/s] of reaction `j` for the counts `n`.
    pub fn propensity(&self, j: usize, n: &[i64]) -> f64 {
        let r = &self.reactions[j];
        for &(s, c) in &r.need {
            if n[s] < c {
                return 0.0;
            }
        }
        match &r.law {
            Law::MassAction { c, reactants } => {
                let mut a = *c;
                for &(s, mult) in reactants {
                    let mut x = n[s];
                    for _ in 0..mult {
                        if x <= 0 {
                            return 0.0;
                        }
                        a *= x as f64;
                        x -= 1;
                    }
                }
                a
            }
            Law::Michaelis {
                k,
                km,
                enzyme,
                substrate,
            } => {
                let cc = |s: &Slot| match s {
                    Slot::Dyn(i) => n[*i].max(0) as f64 / self.nav,
                    Slot::Fixed(c) => *c,
                };
                let (e, s) = (cc(enzyme), cc(substrate));
                if e <= 0.0 || s <= 0.0 {
                    return 0.0;
                }
                k * e * s / (km + s) * self.nav
            }
        }
    }

    /// Molecules of one unit of state [µmol/L].
    pub fn molecules_per_umol(&self) -> f64 {
        1e-6 * self.nav
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;

    pub(crate) fn mat_of(toml: &str) -> (SimEnv, SimMatrices) {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        (env, mat)
    }

    const SELF: &str = r#"
        [[reactions.k_reaction]]
        reaction = "2 OH_r -> H2O2"
        k_value = 5.5e9
        [[reactions.k_reaction]]
        reaction = "H_r + H_r -> H2"
        k_value = 7.8e9
        [[reactions.k_reaction]]
        reaction = "e_aq + OH_r -> OH_minus"
        k_value = 3.0e10
        [[reactions.k_reaction]]
        reaction = "e_aq + H2O -> H_r + OH_minus"
        k_value = 19.0
    "#;

    #[test]
    fn self_reaction_uses_n_times_n_minus_one() {
        let (env, mat) = mat_of(SELF);
        let v = 1e-15;
        let m = SsaModel::new(&mat, v).unwrap();
        let nav = AVOGADRO * v;
        let idx = |l: &str| env.species.index_of(l).unwrap();
        let mut n = vec![0i64; m.n_dyn];
        n[idx("OH_r")] = 10;
        n[idx("H_r")] = 7;
        n[idx("e_aq")] = 3;
        // 2 OH_r -> H2O2 and H_r + H_r -> H2: k / (N_A V) * n (n - 1), with the
        // k of the rate law r = k [A]^2 (a scalar k_value of a single-reactant
        // reaction is divided by its stoichiometry when the config is read).
        assert!((mat.k.k_values[0] - 5.5e9 / 2.0).abs() <= 1e-9 * 5.5e9);
        let a0 = m.propensity(0, &n);
        assert!((a0 - mat.k.k_values[0] / nav * 90.0).abs() <= 1e-12 * a0);
        let a1 = m.propensity(1, &n);
        assert!((a1 - mat.k.k_values[1] / nav * 42.0).abs() <= 1e-12 * a1);
        for j in [0, 1] {
            let sp = if j == 0 { idx("OH_r") } else { idx("H_r") };
            let d = m.reactions[j].delta.iter().find(|(s, _)| *s == sp).unwrap();
            assert_eq!(d.1, -2, "two molecules removed per firing");
            assert_eq!(m.reactions[j].order, 2);
        }
        // A + B: k / (N_A V) * nA * nB.
        let a2 = m.propensity(2, &n);
        assert!((a2 - 3.0e10 / nav * 30.0).abs() <= 1e-12 * a2);
        // Fixed reactant: k * [H2O] * n, first order in the dynamic species.
        let water = mat.cst_sp[idx("H2O") - mat.n_dyn];
        assert!(water > 50.0);
        let a3 = m.propensity(3, &n);
        assert!((a3 - 19.0 * water * 3.0).abs() <= 1e-12 * a3);
        assert_eq!(m.reactions[3].order, 1);
        // One molecule cannot react with itself.
        n[idx("OH_r")] = 1;
        assert_eq!(m.propensity(0, &n), 0.0);
    }

    #[test]
    fn mean_field_limit_matches_the_rate_law() {
        // a / (N_A V) -> k [A]^2 as the count grows.
        let (env, mat) = mat_of(SELF);
        let v = 1e-9;
        let m = SsaModel::new(&mat, v).unwrap();
        let mut n = vec![0i64; m.n_dyn];
        let cc = 2e-6;
        n[env.species.index_of("OH_r").unwrap()] = (cc * m.nav).round() as i64;
        let rate = m.propensity(0, &n) / m.nav;
        assert!((rate - mat.k.k_values[0] * cc * cc).abs() <= 1e-6 * rate);
    }

    #[test]
    fn negative_product_is_a_hard_condition() {
        let (env, mat) = mat_of(
            r#"
            [[reactions.k_reaction]]
            reaction = "OH_r -> H2O2 - O2"
            k_value = 2.0
            "#,
        );
        let m = SsaModel::new(&mat, 1e-15).unwrap();
        let idx = |l: &str| env.species.index_of(l).unwrap();
        let mut n = vec![0i64; m.n_dyn];
        n[idx("OH_r")] = 4;
        assert_eq!(m.propensity(0, &n), 0.0, "no O2 to consume");
        n[idx("O2")] = 1;
        assert_eq!(m.propensity(0, &n), 8.0);
        assert!(m.readers[idx("O2")].contains(&0));
    }

    #[test]
    fn affects_lists_cover_the_readers() {
        let (env, mat) = mat_of(SELF);
        let m = SsaModel::new(&mat, 1e-15).unwrap();
        let oh = env.species.index_of("OH_r").unwrap();
        // Reaction 2 consumes e_aq and OH_r: reactions 0, 2 and 3 read them.
        assert_eq!(m.affects[2], vec![0, 2, 3]);
        assert_eq!(m.readers[oh], vec![0, 2]);
    }

    #[test]
    fn reversible_pairs_and_their_bimolecular_direction() {
        // The radiopyo network: one direction of every pair (declared or
        // inferred) is bimolecular once H2O (fixed) is dropped.
        let cfg = Config::load("data/CONFIG_PAPER_01.toml").unwrap();
        let env = SimEnv::from_config(&cfg).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let m = SsaModel::new(&mat, 1e-18).unwrap();
        let text = env.reaction_strings();
        assert_eq!(m.pairs, mat.reversible_pairs);
        let bimolecular: Vec<&str> = m
            .pairs
            .iter()
            .map(|p| {
                let (f, b) = (p.forward, p.backward);
                assert_ne!(m.is_bimolecular(f), m.is_bimolecular(b), "{}", text[f]);
                text[if m.is_bimolecular(f) { f } else { b }].as_str()
            })
            .collect();
        assert_eq!(
            bimolecular,
            [
                "H3O_plus + HO2_minus <-> H2O2",
                "O_minus + H3O_plus <-> OH_r",
                "O2_minus + H3O_plus <-> HO2",
                "HO2 + OH_minus <-> O2_minus + H2O",
                "H3O_plus + OH_minus <-> H2O",
                "H2O2 + OH_minus <-> HO2_minus + H2O",
                "H_r + OH_minus <-> e_aq + H2O",
                "e_aq + H3O_plus <-> H_r",
                "OH_r + OH_minus <-> O_minus + H2O",
                // Inferred from one-way reactions.
                "OH_r + H2 -> H_r + H2O",
                "O_minus + O2 -> O3_minus",
            ]
        );
    }

    #[test]
    fn rejects_bad_volume() {
        let (_, mat) = mat_of(SELF);
        assert!(SsaModel::new(&mat, 0.0).is_err());
        assert!(SsaModel::new(&mat, f64::NAN).is_err());
    }
}
