"""Compare radkin's IRT stage with Geant4-DNA IRT_syn on the same events.

Input: a Geant4 reference directory with
  1ps/molecules.csv   the molecules of every event when the chemistry starts
                      (species,x,y,z,event_id,time_ps; positions in nm),
  gt/gt.csv           the molecule counts of the G(t) pass (event_id,species,time_s,count),
  gt/edep.csv         the deposited energy (event_id,edep_eV).
Each event is one radkin run (`radkin irt`) of its own molecules, with several
trajectories, recorded at Geant4's checkpoint times.

For a species and a time, N is the number of molecules summed over the events and
G = N / (sum E / 100 eV). The two codes use the same tracks but independent random
numbers, so the difference is tested against two noises:

  z = (N_radkin - N_geant4) / sqrt(sum_e s_e^2 (1 + 1/T) + w^2 sum_e d_e^2)

- the chemistry: s_e is the standard deviation of radkin's T trajectories of event e
  (taken as the noise of Geant4's single realisation too);
- the initial state: the seed pass and the G(t) pass of Geant4 are two realisations of
  the pre-chemical stage of the same track. d_e is the difference of their molecule
  counts at the start, and w = min(1, N_geant4(t) / N_geant4(start)) the fraction of
  the species still there.
The variance has a floor of one molecule.

  python -I tools/compare_irt_geant4.py --events 100 --boundary outflow

The chemistry is data/chem_g4.toml. The report is printed and written to OUT/report.md.
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G4_DEFAULT = os.path.join(ROOT, ".scratch", "tests",
                          "2026-11-15__geant4-single-track-gvalues-10us", "g4")
OUT_DEFAULT = os.path.join(ROOT, ".scratch", "tests", "2026-11-15__irt-syn-port",
                           "compare_irt_geant4")

# Geant4 user IDs -> radkin labels (src/voxel/seed.rs DEFAULT_GEANT4_ALIASES)
NAMES = {
    "e_aq": "e_aq", "°OH": "OH_r", "H3Op": "H3O_plus", "H": "H_r", "H2": "H2",
    "H2O2": "H2O2", "OHm": "OH_minus", "HO2°": "HO2", "HO2m": "HO2_minus", "O2": "O2",
    "O2m": "O2_minus", "O3": "O3", "O3m": "O3_minus", "Om": "O_minus", "Oxy": "O_atom",
}
MAIN = ["e_aq", "OH_r", "H_r", "H3O_plus", "OH_minus", "H2O2", "H2"]
DECADES = [1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]


def read_geant4(g4):
    edep = {}
    with open(os.path.join(g4, "gt", "edep.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            edep[int(r["event_id"])] = float(r["edep_eV"])
    counts = {}  # (event, label, time) -> count
    times = set()
    with open(os.path.join(g4, "gt", "gt.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            t = float(r["time_s"])
            times.add(t)
            counts[(int(r["event_id"]), NAMES[r["species"]], t)] = int(r["count"])
    t_start = None
    with open(os.path.join(g4, "1ps", "molecules.csv"), encoding="utf-8") as f:
        t_start = float(next(csv.DictReader(f))["time_ps"]) * 1e-12
    return edep, counts, sorted(times), t_start


def write_config(path, args, event, energy, times, t_start, t_end):
    mol = os.path.join(args.g4, "1ps", "molecules.csv").replace("\\", "/")
    chem = os.path.join(ROOT, "data", "chem_g4.toml").replace("\\", "/")
    lines = [
        "[includes]", f'toml = ["{chem}"]', "",
        "[irt]", f"t_end = {t_end!r}", f"trajectories = {args.trajectories}",
        f"random_seed = {args.seed + event}", f"start_contact = {str(args.start_contact).lower()}",
        f'boundary = "{args.boundary}"',
        "record_times = [" + ", ".join(repr(t) for t in times) + "]",
    ] + args.irt_key + [
        "", "[irt.box]", f"lo_nm = {-args.half_box}", f"hi_nm = {args.half_box}",
        "", "[[irt.seed]]", f'file = "{mol}"', f"time = {t_start!r}", f"events = [{event}]",
        f"energy_ev = {energy!r}", "",
    ]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))


def run_event(args, event, energy, times, t_start, t_end):
    d = os.path.join(args.out, "runs", f"e{event:03d}")
    os.makedirs(d, exist_ok=True)
    cfg = os.path.join(d, "config.toml")
    write_config(cfg, args, event, energy, times, t_start, t_end)
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    p = subprocess.run([args.radkin, "irt", cfg, "--out", os.path.join(d, "out")],
                       capture_output=True, text=True, env=env)
    if p.returncode != 0:
        raise SystemExit(f"radkin irt failed on event {event}:\n{p.stderr}")
    with open(os.path.join(d, "out", "yields.csv"), encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    with open(os.path.join(d, "out", "irt.json"), encoding="utf-8") as f:
        meta = json.load(f)
    return event, rows, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g4", default=G4_DEFAULT, help="Geant4 reference directory")
    ap.add_argument("--out", default=OUT_DEFAULT, help="output directory")
    ap.add_argument("--radkin", default=os.path.join(ROOT, "target", "release", "radkin.exe"
                                                     if os.name == "nt" else "radkin"))
    ap.add_argument("--events", type=int, default=100, help="number of events (lowest ids)")
    ap.add_argument("--trajectories", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--boundary", choices=["none", "outflow"], default="none")
    ap.add_argument("--half-box", type=float, default=1000.0,
                    help="half edge of the Geant4 water box [nm]")
    ap.add_argument("--start-contact", action="store_true",
                    help="apply the contact rule at the start (off: Geant4 already did at 1 ps)")
    ap.add_argument("--irt-key", action="append", default=[],
                    help='extra line of the [irt] table, e.g. --irt-key "cutoff_probability = 0.001"')
    ap.add_argument("--t-max", type=float, default=None,
                    help="compare up to this time only [s]")
    ap.add_argument("--limit", type=float, default=3.0, help="allowed |z|")
    ap.add_argument("--trace", action="append", default=[],
                    help="also print every compared time of this species (radkin label)")
    args = ap.parse_args()

    edep, g4, times, t_start = read_geant4(args.g4)
    events = sorted(edep)[: args.events]
    t_end = max(times)
    rec = [t for t in times if t_start < t < t_end]
    os.makedirs(args.out, exist_ok=True)
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        runs = list(pool.map(lambda e: run_event(args, e, edep[e], rec, t_start, t_end), events))

    energy = sum(edep[e] for e in events) / 100.0
    species = [k[2:] for k in runs[0][1][0] if k.startswith("N_") and not k.endswith("_std")]
    T = args.trajectories
    cmp_times = [t for t in rec + [t_end] if args.t_max is None or t <= args.t_max]
    table = {}  # (species, t) -> (N_radkin, N_g4, sd)
    for _, rows, _ in runs:
        by_t = {float(r["t"]): r for r in rows}
        missing = [t for t in cmp_times if t not in by_t]
        if missing:
            raise SystemExit(f"record times missing from a radkin run: {missing[:3]}")
    t_first = min(times)
    indexed = [(event, {float(r["t"]): r for r in rows}) for event, rows, _ in runs]
    for sp in species:
        # Difference of the two initial states: Geant4's first checkpoint against the
        # seed molecules radkin starts from.
        init_var = n0 = 0.0
        for event, by_t in indexed:
            g0 = g4.get((event, sp, t_first), 0)
            init_var += (g0 - float(by_t[t_start]["N_" + sp])) ** 2
            n0 += g0
        for t in cmp_times:
            n_rk = var = n_g4 = 0.0
            for event, by_t in indexed:
                r = by_t[t]
                n_rk += float(r["N_" + sp])
                var += float(r["N_" + sp + "_std"]) ** 2 * (1.0 + 1.0 / T)
                n_g4 += g4.get((event, sp, t), 0)
            survive = min(1.0, n_g4 / n0) if n0 > 0 else 1.0
            var += init_var * survive ** 2
            table[(sp, t)] = (n_rk, n_g4, math.sqrt(max(var, 1.0)))

    def z(sp, t):
        n_rk, n_g4, sd = table[(sp, t)]
        return (n_rk - n_g4) / sd

    absorbed = {}
    reactions = {}
    for _, _, meta in runs:
        for k, v in meta["absorbed"].items():
            absorbed[k] = absorbed.get(k, 0) + v / T
        for k, v in meta["reactions"].items():
            reactions[k] = reactions.get(k, 0) + v / T
    seeded = sum(m["seed_molecules"] for _, _, m in runs)

    out = []
    out.append(f"# radkin IRT stage against Geant4-DNA IRT_syn: {len(events)} events\n")
    out.append(f"- Geant4 reference: `{os.path.relpath(args.g4, ROOT)}`; chemistry `data/chem_g4.toml`.")
    out.append(f"- radkin: full resampling, {T} trajectories per event, boundary `{args.boundary}`"
               + (f" (box of half edge {args.half_box:g} nm)" if args.boundary != "none" else "")
               + f", start_contact = {str(args.start_contact).lower()}"
               + ("".join(f", `{k}`" for k in args.irt_key)) + ".")
    out.append(f"- {seeded} molecules at {t_start:.4g} s, {energy * 100:.1f} eV deposited; "
               f"compared at {len(cmp_times)} times up to {max(cmp_times):.3g} s.")
    out.append(f"- Reactions per run of all events: {sum(reactions.values()):.1f}; "
               f"molecules that left the box: {sum(absorbed.values()):.1f}.\n")
    out.append("G in molecules / 100 eV: radkin / Geant4 (z).\n")
    dec = [min(cmp_times, key=lambda t: abs(math.log(t / d))) for d in DECADES
           if args.t_max is None or d <= args.t_max * 1.001]
    out.append("| species | " + " | ".join(f"{t:.0e} s" for t in dec) + " |")
    out.append("|---|" + "---|" * len(dec))
    shown = [s for s in MAIN if s in species]
    shown += [s for s in species if s not in shown
              and max(table[(s, t)][1] for t in cmp_times) >= 10]
    for sp in shown:
        cells = []
        for t in dec:
            n_rk, n_g4, _ = table[(sp, t)]
            cells.append(f"{n_rk / energy:.3f} / {n_g4 / energy:.3f} ({z(sp, t):+.1f})")
        out.append(f"| {sp} | " + " | ".join(cells) + " |")
    out.append("")
    worst = sorted(((abs(z(sp, t)), sp, t) for sp in species for t in cmp_times), reverse=True)
    n_tests = len(species) * len(cmp_times)
    over = [w for w in worst if w[0] > args.limit]
    out.append(f"Largest |z| over {len(species)} species x {len(cmp_times)} times ({n_tests} values):")
    for a, sp, t in worst[:8]:
        n_rk, n_g4, sd = table[(sp, t)]
        out.append(f"- {sp} at {t:.3e} s: z = {z(sp, t):+.2f} (radkin {n_rk:.1f}, Geant4 {n_g4:.0f}, sd {sd:.1f})")
    rms = math.sqrt(sum(w[0] ** 2 for w in worst) / n_tests)
    out.append(f"\nRMS of z: {rms:.2f} (1 expected). Values above {args.limit:g}: {len(over)} of {n_tests}.")
    dec_over = [(sp, t) for sp in shown for t in dec if abs(z(sp, t)) > args.limit]
    ok = not dec_over
    out.append(f"\nAt the decade times, species of the table: "
               + ("every value within " if ok else f"{len(dec_over)} values outside ")
               + f"{args.limit:g} standard errors"
               + ("." if ok else ": " + ", ".join(f"{sp} at {t:.0e} s" for sp, t in dec_over) + "."))
    if absorbed:
        out.append("\nMolecules that left the box, per run of all events: "
                   + ", ".join(f"{k} {v:.1f}" for k, v in sorted(absorbed.items())) + ".")
    for sp in args.trace:
        out.append(f"\nEvery compared time of {sp}: t [s], N radkin, N Geant4, z.\n")
        out.append("```")
        for t in cmp_times:
            n_rk, n_g4, _ = table[(sp, t)]
            out.append(f"{t:.4e} {n_rk:9.1f} {n_g4:7.0f} {z(sp, t):+6.2f}")
        out.append("```")
    text = "\n".join(out) + "\n"
    with open(os.path.join(args.out, "report.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    sys.stdout.reconfigure(encoding="utf-8")
    print(text)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
