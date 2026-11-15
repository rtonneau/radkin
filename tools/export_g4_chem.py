"""Turn a Geant4-DNA chemistry dump into a radkin chemistry TOML.

The dump is the reaction table Geant4 used at run time (species, D, charge, radius,
k_obs, type, radii, products), written by TimeStepAction::DumpChemistry of
${DEV_DIR}/GEANT4/SIM/radkin. The one in the repo is from Geant4-DNA 11.4.3,
chemistry option3.

Mapping rules:
- Second-order reactions: k_obs in M-1 s-1, and the Geant4 reaction type as
  `irt_type` ("full" for type 0, "partial" for type 1), which the IRT stage reads.
- Self-reactions "A + A": Geant4's k_obs is the rate constant of the reaction events
  (events = k_obs [A]^2). radkin divides a scalar k_value of "2 A -> ..." by 2 (the
  n.k convention), so the TOML holds 2 k_obs.
- "X + NoneM": first order, k in s-1.
- "X + H2O(B)": pseudo-first order with water: k = k1 / [H2O] (radkin multiplies
  the constant species H2O into the rate).
- "X + H3Op(B)" / "X + OHm(B)": the pH 7 bulk ions. In radkin bulk and track ions are
  one species, so the matching track reaction is kept and the bulk comes from pH = 7.
  Each one is checked: a track reaction with the same reactants and products must exist.
- Water autoionisation is added to keep pH 7 against H3O+ + OH- (Geant4's bulk is
  constant; radkin would otherwise neutralise it).
- Deaerated: no O2 initially.
Diffusion coefficients go to a JSON next to the TOML (a [voxel.diffusion] table).

  python -I tools/export_g4_chem.py [DUMP [OUT.toml]]

Defaults: tests/data/irt/g4_chemistry_dump.txt and data/chem_g4.toml.
"""

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H2O = 55.93  # mol/L, as in data/voxel_example.toml

# Geant4 user IDs -> radkin labels (src/voxel/seed.rs DEFAULT_GEANT4_ALIASES)
NAMES = {
    "e_aq": "e_aq", "°OH": "OH_r", "H3Op": "H3O_plus", "H": "H_r", "H2": "H2",
    "H2O2": "H2O2", "OHm": "OH_minus", "HO2°": "HO2", "HO2m": "HO2_minus", "O2": "O2",
    "O2m": "O2_minus", "O3": "O3", "O3m": "O3_minus", "Om": "O_minus", "Oxy": "O_atom",
}
BULK_IONS = {"H3Op(B)": "H3Op", "OHm(B)": "OHm"}
IRT_TYPES = {0: "full", 1: "partial"}


def parse(path):
    species, reactions = {}, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "species":
                species[p[2]] = dict(name=p[1], D=float(p[3]), charge=int(p[4]))
            elif p[0] == "reaction":
                reactions.append(dict(id=int(p[1]), a=p[2], b=p[3], type=int(p[4]),
                                      k2=float(p[6]), k1=float(p[7]),
                                      products=[] if p[10] == "-" else p[10].split(";")))
    return species, reactions


def side(names):
    """radkin side of a reaction, with stoichiometric coefficients."""
    out, seen = [], {}
    for n in names:
        seen[n] = seen.get(n, 0) + 1
    for n, c in seen.items():
        out.append(f"{c} {n}" if c > 1 else n)
    return " + ".join(out) if out else "H2O"


def label(uid):
    if uid not in NAMES:
        raise SystemExit(f"no radkin label for Geant4 species {uid!r}")
    return NAMES[uid]


def key(r):
    return (frozenset([r["a"], r["b"]]), tuple(sorted(r["products"])))


def main(dump, out):
    species, reactions = parse(dump)
    track = {key(r) for r in reactions}
    lines = [
        "# Geant4-DNA chemistry (11.4.3, option3) exported from the reaction table Geant4",
        f"# used at run time: {os.path.basename(dump)} by tools/export_g4_chem.py.",
        "# Deaerated, pH 7 background. Rates are the observed rate constants k_obs;",
        "# a self-reaction holds 2 k_obs (radkin divides it by 2).",
        "# irt_type is the Geant4 reaction type: \"full\" (0) or \"partial\" (1).",
        "",
        "[parameters]",
        "pH = 7",
        "force_cst_pH = false",
        "",
        "[concentrations.initial]",
        "",
    ]
    dropped = []
    for r in reactions:
        a, b = r["a"], r["b"]
        prods = [label(p) for p in r["products"]]
        irt = None
        if b in BULK_IONS or a in BULK_IONS:
            ion, other = (b, a) if b in BULK_IONS else (a, b)
            twin = (frozenset([other, BULK_IONS[ion]]), tuple(sorted(r["products"])))
            if twin not in track:
                raise SystemExit(f"bulk reaction {r['id']} ({a} + {b}) has no track twin")
            dropped.append(r["id"])
            continue
        if b == "NoneM" or a == "NoneM":
            x = a if b == "NoneM" else b
            rx, kv, note = f"{label(x)} -> {side(prods)}", r["k1"], "first order, s-1"
        elif b == "H2O(B)" or a == "H2O(B)":
            x = a if b == "H2O(B)" else b
            rx = f"{label(x)} + H2O -> {side(prods)}"
            kv, note = r["k1"] / H2O, f"{r['k1']:.4g} s-1 / [H2O] {H2O}"
        else:
            rx = f"{side([label(a), label(b)])} -> {side(prods)}"
            irt = IRT_TYPES[r["type"]]
            if a == b:
                kv, note = 2 * r["k2"], f"M-1 s-1, 2 x k_obs {r['k2']:.6g}"
            else:
                kv, note = r["k2"], "M-1 s-1"
        # "2 O_atom" would parse as a species: write the self-reaction of O_atom long-hand
        rx = rx.replace("2 O_atom ->", "O_atom + O_atom ->")
        lines += ["[[reactions.k_reaction]]", f'reaction = "{rx}"',
                  f"k_value = {kv:.6g}   # G4 reaction {r['id']} ({note})"]
        if irt:
            lines.append(f'irt_type = "{irt}"')
        lines.append("")
    k_auto = 1.13e11 * 1e-14 / H2O
    lines += ["# not in Geant4 (its bulk ions are constant): keeps pH 7 against H3O+ + OH-",
              "[[reactions.k_reaction]]", 'reaction = "H2O -> H3O_plus + OH_minus"',
              f"k_value = {k_auto:.6g}", "",
              "# Bulk (B) ion reactions dropped (track twins kept): "
              + ", ".join(map(str, dropped))]
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    diff = {NAMES[u]: s["D"] for u, s in species.items() if u in NAMES}
    with open(os.path.splitext(out)[0] + ".diffusion.json", "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(diff, f, indent=1)
        f.write("\n")
    print(f"{out}: {sum(1 for l in lines if l.startswith('reaction ='))} reactions, "
          f"{len(dropped)} bulk dropped, D for {len(diff)} species")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0] if args else os.path.join(ROOT, "tests", "data", "irt", "g4_chemistry_dump.txt"),
         args[1] if len(args) > 1 else os.path.join(ROOT, "data", "chem_g4.toml"))
