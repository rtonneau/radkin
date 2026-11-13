//! Reaction-string parser: a hand-written port of radiopyo's Lark grammar
//! (`parser/grammar/reaction_grammar.lark`) and its transformer.
//!
//! - `A + B -> C`: k-reaction; `-> C`: constant production;
//!   `A <-> B + C`: equilibrium; `A -- enzyme >> B`: enzymatic.
//! - A species is `[n]formula` or a plain name (`catalase`, `trash`).
//!   Formulas are atoms `Capital[lower][digits]`, then `_r`, then charges
//!   (`_plus`, `_minus`, `+`, `-`); `e_aq` is the hydrated electron.
//! - Reactants with the same label are summed; products are summed with
//!   their sign, a net negative stoichiometry is a [`Operation::Subtract`]
//!   ("negative" product), a net zero one is dropped.

use std::collections::BTreeSet;
use std::fmt;

use indexmap::IndexMap;

use crate::error::{RadkinError, Result};

/// Sign of a species in a reaction side.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Operation {
    Add,
    Subtract,
}

impl Operation {
    /// Name used in radiopyo's reaction keys (`ope.name`).
    pub fn name(&self) -> &'static str {
        match self {
            Self::Add => "ADD",
            Self::Subtract => "SUBTRACT",
        }
    }
}

impl fmt::Display for Operation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Add => "+",
            Self::Subtract => "-",
        })
    }
}

/// Kind of reaction, from its arrow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionKind {
    K,
    Constant,
    Equilibrium,
    Enzymatic,
}

/// Species side: label -> (stoichiometry, operation), in order of appearance.
pub type Side = IndexMap<String, (u32, Operation)>;

/// A parsed reaction string (radiopyo `GenericReaction` without constants).
#[derive(Debug, Clone, PartialEq)]
pub struct GenericReaction {
    pub kind: ReactionKind,
    pub reactants: Side,
    pub products: Side,
    pub enzyme: Option<String>,
}

/// Order-insensitive identity of a reaction (radiopyo `reaction_key`).
pub type ReactionKey = (
    BTreeSet<(String, u32)>,
    BTreeSet<(String, u32, &'static str)>,
    Option<String>,
);

fn err(input: &str, msg: impl Into<String>) -> RadkinError {
    RadkinError::Parse {
        input: input.to_string(),
        msg: msg.into(),
    }
}

/// Parse a reaction string.
pub fn parse_reaction(s: &str) -> Result<GenericReaction> {
    let text = s.trim();
    let (kind, left, right, enzyme) = if let Some(i) = text.find("<->") {
        (ReactionKind::Equilibrium, &text[..i], &text[i + 3..], None)
    } else if let (Some(i), Some(j)) = (text.find("--"), text.find(">>")) {
        if j < i {
            return Err(err(s, "'>>' before '--'"));
        }
        let enzyme = text[i + 2..j].trim();
        if enzyme.is_empty() || !enzyme.chars().all(is_word) {
            return Err(err(s, format!("invalid enzyme name '{enzyme}'")));
        }
        (
            ReactionKind::Enzymatic,
            &text[..i],
            &text[j + 2..],
            Some(enzyme.to_string()),
        )
    } else if let Some(i) = text.find("->") {
        let kind = if text[..i].trim().is_empty() {
            ReactionKind::Constant
        } else {
            ReactionKind::K
        };
        (kind, &text[..i], &text[i + 2..], None)
    } else {
        return Err(err(s, "no reaction arrow ('->', '<->' or '-- enzyme >>')"));
    };
    let left_items = parse_side(s, left, false)?;
    let right_items = parse_side(s, right, true)?;
    if kind == ReactionKind::Enzymatic && left_items.len() != 1 {
        return Err(err(s, "an enzymatic reaction has exactly one substrate"));
    }
    if matches!(
        kind,
        ReactionKind::K | ReactionKind::Equilibrium | ReactionKind::Enzymatic
    ) && left_items.is_empty()
    {
        return Err(err(s, "no reactant"));
    }
    let mut reactants = Side::new();
    for (label, stoi, _) in left_items {
        reactants.entry(label).or_insert((0, Operation::Add)).0 += stoi;
    }
    let mut sums: IndexMap<String, i64> = IndexMap::new();
    for (label, stoi, ope) in right_items {
        let signed = if ope == Operation::Add {
            stoi as i64
        } else {
            -(stoi as i64)
        };
        *sums.entry(label).or_insert(0) += signed;
    }
    let products = sums
        .into_iter()
        .filter(|(_, v)| *v != 0)
        .map(|(k, v)| {
            let ope = if v < 0 {
                Operation::Subtract
            } else {
                Operation::Add
            };
            (k, (v.unsigned_abs() as u32, ope))
        })
        .collect();
    Ok(GenericReaction {
        kind,
        reactants,
        products,
        enzyme,
    })
}

/// Parse a single species name and return its label.
pub fn parse_species(s: &str) -> Result<String> {
    let text = s.trim();
    let items = parse_side(s, text, false)?;
    match items.as_slice() {
        [(label, 1, Operation::Add)] => Ok(label.clone()),
        _ => Err(err(s, "not a single chemical species")),
    }
}

/// Swap reactants and products (radiopyo `revert_reaction`).
pub fn revert(g: &GenericReaction) -> Result<GenericReaction> {
    if g.products.values().any(|(_, o)| *o != Operation::Add) {
        return Err(RadkinError::Value(
            "Cannot revert reaction with subtraction operation".to_string(),
        ));
    }
    Ok(GenericReaction {
        kind: g.kind,
        reactants: g.products.clone(),
        products: g.reactants.clone(),
        enzyme: g.enzyme.clone(),
    })
}

/// Order-insensitive key of a reaction string (radiopyo `reaction_key`).
pub fn reaction_key(s: &str) -> Result<ReactionKey> {
    let g = parse_reaction(s)?;
    Ok((
        g.reactants
            .iter()
            .map(|(k, (n, _))| (k.clone(), *n))
            .collect(),
        g.products
            .iter()
            .map(|(k, (n, o))| (k.clone(), *n, o.name()))
            .collect(),
        g.enzyme,
    ))
}

fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// Parse one side into (label, stoichiometry, sign) items, in order.
/// `signed`: a leading sign is allowed (products) and `-` separates.
fn parse_side(input: &str, side: &str, signed: bool) -> Result<Vec<(String, u32, Operation)>> {
    let chars: Vec<char> = side.chars().collect();
    let mut i = 0;
    let mut out = Vec::new();
    let mut pending_sign: Option<Operation> = None;
    let skip_ws = |i: &mut usize| {
        while *i < chars.len() && chars[*i].is_whitespace() {
            *i += 1;
        }
    };
    loop {
        skip_ws(&mut i);
        if i >= chars.len() {
            if pending_sign.is_some() {
                return Err(err(input, "a sign is not followed by a species"));
            }
            break;
        }
        let c = chars[i];
        if c == '+' || c == '-' {
            if pending_sign.is_some() {
                return Err(err(input, "two signs in a row"));
            }
            if out.is_empty() && !signed {
                return Err(err(input, format!("unexpected '{c}'")));
            }
            if c == '-' && !signed {
                return Err(err(input, "'-' is only allowed between products"));
            }
            pending_sign = Some(if c == '+' {
                Operation::Add
            } else {
                Operation::Subtract
            });
            i += 1;
            continue;
        }
        if !is_word(c) {
            return Err(err(input, format!("unexpected character '{c}'")));
        }
        if !out.is_empty() && pending_sign.is_none() {
            return Err(err(input, "species not separated by '+'"));
        }
        let start = i;
        while i < chars.len() && is_word(chars[i]) {
            i += 1;
        }
        let mut word: String = chars[start..i].iter().collect();
        if word.chars().all(|c| c.is_ascii_digit()) {
            // Multiplicity separated from its species by blanks ("2 H2O").
            skip_ws(&mut i);
            let s2 = i;
            while i < chars.len() && is_word(chars[i]) {
                i += 1;
            }
            if s2 == i {
                return Err(err(
                    input,
                    format!("multiplicity '{word}' without a species"),
                ));
            }
            word.extend(&chars[s2..i]);
        }
        // Attached +/- run: charges, except a last char that separates the next species.
        let run_start = i;
        while i < chars.len() && (chars[i] == '+' || chars[i] == '-') {
            i += 1;
        }
        let mut run: Vec<char> = chars[run_start..i].to_vec();
        // Next non-blank char starts a species ("H_r+ OH_minus", "A+B").
        let mut k = i;
        while k < chars.len() && chars[k].is_whitespace() {
            k += 1;
        }
        let next_starts_species = k < chars.len() && is_word(chars[k]);
        let mut separator = None;
        if !run.is_empty() && next_starts_species {
            separator = run.pop();
        }
        let (label, stoi) = species_label(input, &word, &run)?;
        out.push((label, stoi, pending_sign.take().unwrap_or(Operation::Add)));
        if let Some(sep) = separator {
            if sep == '-' && !signed {
                return Err(err(input, "'-' is only allowed between products"));
            }
            pending_sign = Some(if sep == '+' {
                Operation::Add
            } else {
                Operation::Subtract
            });
        }
    }
    Ok(out)
}

/// Label and multiplicity of a species token (`word` + attached charge chars).
fn species_label(input: &str, word: &str, charges: &[char]) -> Result<(String, u32)> {
    let digits: String = word.chars().take_while(|c| c.is_ascii_digit()).collect();
    if !digits.is_empty() && digits.len() < word.len() {
        let rest = &word[digits.len()..];
        if let Some(label) = formula_label(rest, charges) {
            let n: u32 = digits.parse().map_err(|_| err(input, "bad multiplicity"))?;
            return Ok((label, n));
        }
    }
    if let Some(label) = formula_label(word, charges) {
        return Ok((label, 1));
    }
    if !charges.is_empty() {
        return Err(err(
            input,
            format!("charge on a non-formula species '{word}'"),
        ));
    }
    Ok((word.to_string(), 1))
}

/// Canonical label of a chemical formula, `None` if `word` is not one.
fn formula_label(word: &str, charges: &[char]) -> Option<String> {
    if word == "e_aq" {
        return charges.is_empty().then(|| "e_aq".to_string());
    }
    let b = word.as_bytes();
    let mut i = 0;
    let mut atoms = String::new();
    let mut n_atoms = 0;
    while i < b.len() && b[i].is_ascii_uppercase() {
        let start = i;
        i += 1;
        if i < b.len() && b[i].is_ascii_lowercase() {
            i += 1;
        }
        let name = &word[start..i];
        let d0 = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        let stoi: u64 = if d0 < i { word[d0..i].parse().ok()? } else { 1 };
        atoms.push_str(name);
        if stoi > 1 {
            atoms.push_str(&stoi.to_string());
        }
        n_atoms += 1;
    }
    if n_atoms == 0 {
        return None;
    }
    let mut rest = &word[i..];
    let mut radical = false;
    if let Some(r) = rest.strip_prefix("_r")
        && (r.is_empty() || r.starts_with('_'))
    {
        radical = true;
        rest = r;
    }
    let mut charge: i32 = 0;
    let mut sign_seen: Option<bool> = None;
    let mut push = |plus: bool| -> Option<()> {
        match sign_seen {
            Some(s) if s != plus => None,
            _ => {
                sign_seen = Some(plus);
                charge += if plus { 1 } else { -1 };
                Some(())
            }
        }
    };
    while !rest.is_empty() {
        let (plus, r) = match rest.strip_prefix("_plus") {
            Some(r) => (true, r),
            None => (false, rest.strip_prefix("_minus")?),
        };
        push(plus)?;
        rest = r;
    }
    for c in charges {
        push(*c == '+')?;
    }
    let mut label = atoms;
    if radical {
        label.push_str("_r");
    }
    let suffix = if charge > 0 { "_plus" } else { "_minus" };
    for _ in 0..charge.unsigned_abs() {
        label.push_str(suffix);
    }
    Some(label)
}
