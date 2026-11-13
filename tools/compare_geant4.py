#!/usr/bin/env python
"""Plot the G(t) of a `radkin voxel` run, optionally against Geant4-DNA yields.

    python -I tools/compare_geant4.py OUT_DIR [--geant4 g4.csv] [--species e_aq,OH_r]
        [--map "e_aq^-1=e_aq,°OH^0=OH_r"] [--save plot.png]

`g4.csv` has columns `t,species,G` (time in s, G in #/100 eV), e.g. exported
from a Geant4-DNA chem example's species scorer. Geant4 names are mapped with
`--map` (the standard water species are mapped by default).
"""

import argparse
import csv
import math
import sys
from pathlib import Path

DEFAULT_MAP = {
    "e_aq^-1": "e_aq", "e_aq": "e_aq",
    "°OH^0": "OH_r", "°OH": "OH_r", "OH^0": "OH_r",
    "H3O^1": "H3O_plus", "H3Op": "H3O_plus",
    "H^0": "H_r", "H": "H_r",
    "H_2^0": "H2", "H2": "H2",
    "H2O2^0": "H2O2", "H2O2": "H2O2",
    "OH^-1": "OH_minus", "OHm": "OH_minus",
}


def read_radkin(out_dir: Path):
    with open(out_dir / "yields.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows or not any(k.startswith("G_") for k in rows[0]):
        sys.exit("yields.csv has no G_ columns: give the seeds an energy_ev")
    t = [float(r["t"]) for r in rows]
    g = {
        k[2:]: [float(r[k]) for r in rows]
        for k in rows[0]
        if k.startswith("G_")
    }
    return t, g


def read_geant4(path: Path, mapping):
    out = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sp = mapping.get(r["species"], r["species"])
            out.setdefault(sp, ([], []))
            out[sp][0].append(float(r["t"]))
            out[sp][1].append(float(r["G"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--geant4", type=Path)
    ap.add_argument("--species", default="e_aq,OH_r,H3O_plus,H_r,H2,H2O2")
    ap.add_argument("--map", default="")
    ap.add_argument("--save", type=Path)
    a = ap.parse_args()

    mapping = dict(DEFAULT_MAP)
    for item in filter(None, a.map.split(",")):
        k, v = item.split("=", 1)
        mapping[k.strip()] = v.strip()
    t, g = read_radkin(a.out_dir)
    g4 = read_geant4(a.geant4, mapping) if a.geant4 else {}
    species = [s for s in a.species.split(",") if s in g]

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, sp in enumerate(species):
        color = f"C{i}"
        pts = [(ti, gi) for ti, gi in zip(t, g[sp]) if ti > 0 and math.isfinite(gi)]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, label=f"radkin {sp}")
        if sp in g4:
            ax.plot(*g4[sp], "o", color=color, ms=3, label=f"Geant4-DNA {sp}")
    ax.set_xscale("log")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("G [#/100 eV]")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if a.save:
        fig.savefig(a.save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
