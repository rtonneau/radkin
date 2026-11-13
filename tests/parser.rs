use radkin::chemistry::parser::{Operation, ReactionKind, parse_reaction, parse_species};
use serde_json::Value;

fn golden() -> serde_json::Map<String, Value> {
    let text = std::fs::read_to_string("tests/golden/parser_labels.json").unwrap();
    serde_json::from_str::<Value>(&text)
        .unwrap()
        .as_object()
        .unwrap()
        .clone()
}

fn side(entries: &indexmap::IndexMap<String, (u32, Operation)>) -> Value {
    Value::Array(
        entries
            .iter()
            .map(|(k, (s, o))| serde_json::json!([k, s, o.to_string()]))
            .collect(),
    )
}

#[test]
fn every_golden_reaction_string_gives_radiopyo_labels() {
    for (s, expected) in golden() {
        let g = parse_reaction(&s).unwrap_or_else(|e| panic!("{s}: {e}"));
        assert_eq!(
            side(&g.reactants),
            expected["reactants"],
            "reactants of {s:?}"
        );
        assert_eq!(side(&g.products), expected["products"], "products of {s:?}");
        let enzyme = g.enzyme.map(Value::String).unwrap_or(Value::Null);
        assert_eq!(enzyme, expected["enzyme"], "enzyme of {s:?}");
    }
}

#[test]
fn net_negative_product_is_a_subtraction() {
    let g = parse_reaction("A-> B - 2B + 3 C").unwrap();
    let p: Vec<_> = g.products.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    assert_eq!(
        p,
        vec![("B", (1, Operation::Subtract)), ("C", (3, Operation::Add))]
    );
}

#[test]
fn enzymatic_and_constant_kinds() {
    let g = parse_reaction("2H2O2 -- catalase >> O2 + 2 H2O").unwrap();
    assert_eq!(g.kind, ReactionKind::Enzymatic);
    assert_eq!(g.enzyme.as_deref(), Some("catalase"));
    assert_eq!(g.reactants.get("H2O2"), Some(&(2, Operation::Add)));
    assert_eq!(
        parse_reaction("-> A + B").unwrap().kind,
        ReactionKind::Constant
    );
    assert_eq!(
        parse_reaction("H2O <-> H3O_plus + OH_minus").unwrap().kind,
        ReactionKind::Equilibrium
    );
    assert_eq!(parse_reaction("A + B -> C").unwrap().kind, ReactionKind::K);
}

#[test]
fn empty_products_accepted_but_dangling_plus_rejected() {
    assert!(parse_reaction("A + B -> ").unwrap().products.is_empty());
    let err = parse_reaction("A + -> B").unwrap_err().to_string();
    assert!(err.contains("A + -> B"), "{err}");
}

#[test]
fn species_labels() {
    assert_eq!(parse_species("OH_r").unwrap(), "OH_r");
    assert_eq!(parse_species("Fe++").unwrap(), "Fe_plus_plus");
    assert_eq!(parse_species("e_aq").unwrap(), "e_aq");
    assert_eq!(parse_species(" catalase ").unwrap(), "catalase");
}
