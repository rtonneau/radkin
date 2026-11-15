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

/// Two mass-action reactions that are the two directions of one equilibrium.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionPair {
    /// Index in `k` of the forward direction (the lower index of an
    /// inferred pair).
    pub forward: usize,
    /// Index in `k` of the backward direction.
    pub backward: usize,
    /// Inferred from two one-way reactions rather than declared with `<->`.
    pub inferred: bool,
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
    /// Every pair declared with `<->` (acid/base, then equilibrium, in
    /// declaration order), then the pairs inferred from one-way reactions
    /// (in index order of their forward direction, see
    /// [`infer_one_way_pairs`]).
    pub reversible_pairs: Vec<ReactionPair>,
    /// One-way reactions with more than one candidate reverse, and their
    /// candidates: none of them is paired.
    pub ambiguous_reverses: Vec<(usize, Vec<usize>)>,
    pub n_species: usize,
    pub n_dyn: usize,
}

fn stoichiometry_entry(col: &mut Vec<(usize, f64)>, species: usize, delta: f64) {
    match col.iter_mut().find(|(s, _)| *s == species) {
        Some(e) => e.1 += delta,
        None => col.push((species, delta)),
    }
}

/// Pairs of one-way mass-action reactions that reverse each other.
///
/// `sides[r]` holds the sorted dynamic reactant and product multisets of
/// reaction `r` (fixed species dropped), or `None` if `r` cannot be paired (a
/// direction of a declared pair, a reaction with a negative product, or one
/// without a net change). Reactions `a` and `b` are candidates of each other
/// when `reactants(a) = products(b)` and `products(a) = reactants(b)`. They are
/// paired only when each is the single candidate of the other; every reaction
/// with more than one candidate is returned with its candidates instead.
pub fn infer_one_way_pairs(
    sides: &[Option<(Vec<usize>, Vec<usize>)>],
) -> (Vec<ReactionPair>, Vec<(usize, Vec<usize>)>) {
    let candidates: Vec<Vec<usize>> = sides
        .iter()
        .enumerate()
        .map(|(a, side)| {
            let Some((ra, pa)) = side else {
                return Vec::new();
            };
            sides
                .iter()
                .enumerate()
                .filter(|(b, other)| {
                    *b != a && other.as_ref().is_some_and(|(rb, pb)| rb == pa && pb == ra)
                })
                .map(|(b, _)| b)
                .collect()
        })
        .collect();
    let mut pairs = Vec::new();
    let mut ambiguous = Vec::new();
    for (a, c) in candidates.iter().enumerate() {
        match c.as_slice() {
            [b] if *b > a && candidates[*b] == [a] => pairs.push(ReactionPair {
                forward: a,
                backward: *b,
                inferred: true,
            }),
            [] | [_] => {}
            _ => ambiguous.push((a, c.clone())),
        }
    }
    (pairs, ambiguous)
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
        // `iter_all_k` puts every pair first, forward then backward.
        let n_declared = 2 * env.reactions.iter_pairs().count();
        let mut sides = Vec::new();
        for (j, r) in env.reactions.iter_all_k().enumerate() {
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
            // Candidate of an inferred pair: a one-way reaction with a net
            // change and no negative product.
            let side = if j < n_declared || !negative.is_empty() || col.is_empty() {
                None
            } else {
                let mut reactants: Vec<usize> =
                    slots.iter().copied().filter(|&i| i < n_dyn).collect();
                let mut products = Vec::new();
                for s in r.products.iter().filter(|s| !sp.is_cst(&s.label)) {
                    products.extend(std::iter::repeat_n(idx(&s.label)?, s.stoi as usize));
                }
                reactants.sort_unstable();
                products.sort_unstable();
                Some((reactants, products))
            };
            sides.push(side);
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
        let (inferred, ambiguous_reverses) = infer_one_way_pairs(&sides);
        let reversible_pairs = (0..n_declared / 2)
            .map(|i| ReactionPair {
                forward: 2 * i,
                backward: 2 * i + 1,
                inferred: false,
            })
            .chain(inferred)
            .collect();
        Ok(Self {
            k,
            radiolytic,
            enzymatic: ez,
            cst_sp,
            overshoot,
            reversible_pairs,
            ambiguous_reverses,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;

    #[test]
    fn reversible_pairs_are_the_declared_pairs() {
        let cfg = Config::load("data/CONFIG_PAPER_01.toml").unwrap();
        let env = SimEnv::from_config(&cfg).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let text = env.reaction_strings();
        // The declared forward reaction, and the backward one with the sides
        // swapped, in declaration order (acid/base, then equilibrium).
        let declared = [
            "H2O2 <-> H3O_plus + HO2_minus",
            "OH_r <-> O_minus + H3O_plus",
            "HO2 <-> O2_minus + H3O_plus",
            "HO2 + OH_minus <-> O2_minus + H2O",
            "H2O <-> H3O_plus + OH_minus",
            "H2O2 + OH_minus <-> HO2_minus + H2O",
            "e_aq + H2O <-> H_r + OH_minus",
            "H_r <-> e_aq + H3O_plus",
            "OH_r + OH_minus <-> O_minus + H2O",
        ];
        // The declared pairs come first; the one-way reactions then give two
        // inferred pairs, each direction the single reverse of the other.
        let inferred = [
            ("H_r + H2O -> H2 + OH_r", "OH_r + H2 -> H_r + H2O"),
            ("O_minus + O2 -> O3_minus", "O3_minus -> O_minus + O2"),
        ];
        assert_eq!(mat.reversible_pairs.len(), declared.len() + inferred.len());
        let (first, rest) = mat.reversible_pairs.split_at(declared.len());
        for (p, (f, b)) in rest.iter().zip(inferred) {
            assert!(p.inferred);
            assert_eq!(
                (text[p.forward].as_str(), text[p.backward].as_str()),
                (f, b)
            );
        }
        for (p, want) in first.iter().zip(declared) {
            let (f, b) = (p.forward, p.backward);
            assert!(!p.inferred);
            assert_eq!(text[f], want);
            let (lhs, rhs) = want.split_once(" <-> ").unwrap();
            assert_eq!(text[b], format!("{rhs} <-> {lhs}"));
        }
        // A one-way reaction is not paired with a direction of a declared
        // pair ("e_aq + H2O -> H_r + OH_minus" is also a k_reaction).
        let one_way = text
            .iter()
            .position(|t| t == "e_aq + H2O -> H_r + OH_minus")
            .unwrap();
        assert!(one_way >= 2 * declared.len());
        assert!(
            mat.reversible_pairs
                .iter()
                .all(|p| p.forward != one_way && p.backward != one_way)
        );
        assert!(mat.ambiguous_reverses.is_empty());
    }

    fn build(toml: &str) -> (Vec<String>, SimMatrices) {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        (env.reaction_strings(), SimMatrices::build(&env).unwrap())
    }

    fn k(reaction: &str) -> String {
        format!(
            "[[reactions.k_reaction]]
reaction = \"{reaction}\"
k_value = 1.0e9
"
        )
    }

    #[test]
    fn one_way_reverses_are_inferred_as_pairs() {
        // The Geant4 export writes water autoionisation as two one-way
        // reactions (chem_g4.toml); H2O is fixed and dropped.
        let toml = [
            k("2 H_r -> H2"),
            k("H3O_plus + OH_minus -> H2O"),
            k("e_aq + H_r -> OH_minus + H2"),
            k("H2O -> H3O_plus + OH_minus"),
        ]
        .concat();
        let (text, mat) = build(&toml);
        assert_eq!(
            mat.reversible_pairs,
            [ReactionPair {
                forward: 1,
                backward: 3,
                inferred: true
            }]
        );
        assert_eq!(text[1], "H3O_plus + OH_minus -> H2O");
        assert!(mat.ambiguous_reverses.is_empty());
    }

    #[test]
    fn inferred_pairs_match_dynamic_species_only() {
        // Written with and without H2O, and with a stoichiometry of 2.
        let toml = [
            k("e_aq + H2O -> H_r + OH_minus"),
            k("H_r + OH_minus -> e_aq"),
            k("2 OH_r -> H2O2"),
            k("H2O2 -> 2 OH_r"),
            // Not a reverse: one OH_r only.
            k("H2O2 -> OH_r + H2O"),
        ]
        .concat();
        let (_, mat) = build(&toml);
        let pairs: Vec<(usize, usize)> = mat
            .reversible_pairs
            .iter()
            .map(|p| (p.forward, p.backward))
            .collect();
        assert_eq!(pairs, [(0, 1), (2, 3)]);
    }

    #[test]
    fn a_reaction_with_several_reverses_is_not_paired() {
        let toml = [
            k("H3O_plus + OH_minus -> H2O"),
            k("H2O -> H3O_plus + OH_minus"),
            k("H3O_plus + OH_minus -> 2 H2O"),
        ]
        .concat();
        let (_, mat) = build(&toml);
        assert!(mat.reversible_pairs.is_empty());
        // Only reaction 1 has two candidates; 0 and 2 have one each.
        assert_eq!(mat.ambiguous_reverses, [(1, vec![0, 2])]);
    }

    #[test]
    fn one_way_reactions_are_not_paired_with_a_declared_direction() {
        let toml = format!(
            "[[reactions.equilibrium]]
reaction = \"H2O <-> H3O_plus + OH_minus\"
             k_forward = 2.12e-5
k_backward = 1.18e11
{}{}",
            k("H3O_plus + OH_minus -> H2O"),
            k("H2O -> H3O_plus + OH_minus")
        );
        let (_, mat) = build(&toml);
        // The two one-way reactions pair with each other, not with the
        // declared directions 0 and 1.
        assert_eq!(
            mat.reversible_pairs,
            [
                ReactionPair {
                    forward: 0,
                    backward: 1,
                    inferred: false
                },
                ReactionPair {
                    forward: 2,
                    backward: 3,
                    inferred: true
                },
            ]
        );
        // A single one-way reverse of a declared direction stays unpaired.
        let toml = format!(
            "[[reactions.equilibrium]]
reaction = \"H2O <-> H3O_plus + OH_minus\"
             k_forward = 2.12e-5
k_backward = 1.18e11
{}",
            k("H3O_plus + OH_minus -> H2O")
        );
        let (_, mat) = build(&toml);
        assert_eq!(mat.reversible_pairs.len(), 1);
        assert!(mat.ambiguous_reverses.is_empty());
    }
}
