#!/usr/bin/env python3
"""
Run radiopyo and radkin on the same inputs and compare their outputs.

Both codes read the same wrapper TOML (written in a temporary folder): it
includes the given config and adds the beam, initial concentrations and
parameters given on the command line, so both use radiopyo's merge rules.
Both report at the same times:
  - mode "run": a shared t_eval (log-spaced, --t-eval-n points);
  - chunked modes: N log-spaced points per chunk (--points-on/--points-off),
    the same grid in both codes (radiopyo `thinned_times`).

Example (radiopyo conda env):

    python -I tools/compare_radiopyo.py data/config_Labarbe.toml --mode run \
        --t-span 0 1 --beam-type constant --beam-param dose_rate=1e3 \
        --beam-param max_dose=10 --let 0.201

Exit code 0 when every species, reaction rate and integral agrees within
--rtol (absolute floors: 1e-10 µmol/L for species, 1e-7 x max for rates,
1e-9 x max for integrals), 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("config", type=Path)
    p.add_argument("--mode", choices=("run", "chunked", "single-pulse"), default="run")
    p.add_argument("--t-span", nargs=2, type=float, required=True)
    p.add_argument("--beam-type", choices=("constant", "singlepulse", "pulsed", "macro_pulsed"))
    p.add_argument("--beam-param", action="append", default=[], metavar="KEY=VALUE")
    p.add_argument("--let", type=float, default=0.201)
    p.add_argument("--initial", action="append", default=[], metavar="SPECIES=MOL_PER_L")
    p.add_argument("--variable-gvalues", action="store_true")
    p.add_argument("--t-eval-n", type=int, default=60)
    p.add_argument("--max-step", type=float)
    p.add_argument("--max-step-on", type=float)
    p.add_argument("--max-step-off", type=float)
    p.add_argument("--max-step-trans", type=float)
    p.add_argument("--max-step-final", type=float)
    p.add_argument("--points-on", type=int, default=10)
    p.add_argument("--points-off", type=int, default=20)
    p.add_argument("--production", default="H2O2")
    p.add_argument("--consumption", default="O2")
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--radiopyo", type=Path, default=Path("D:/DEV/PYTHON/radiopyo"))
    p.add_argument("--radkin", type=Path, default=ROOT / "target" / "release" / "radkin")
    p.add_argument("--keep", type=Path, help="keep the radkin outputs in this folder")
    return p.parse_args()


def _kv(items: list[str]) -> dict[str, float]:
    out = {}
    for item in items:
        k, v = item.split("=", 1)
        out[k.strip()] = float(v)
    return out


def write_wrapper(args: argparse.Namespace, folder: Path) -> Path:
    config = args.config.resolve().as_posix()
    lines = ["[includes]", f'toml = ["{config}"]', "", "[parameters]"]
    if args.variable_gvalues:
        lines.append("variable_GValues = true")
    initial = _kv(args.initial)
    if initial:
        lines += ["", "[concentrations.initial]"] + [f"{k} = {v!r}" for k, v in initial.items()]
    if args.beam_type:
        lines += ["", f"[[beam.{args.beam_type}]]"]
        lines += [f"{k} = {v!r}" for k, v in _kv(args.beam_param).items()]
        lines.append(f"LET = {args.let!r}")
    path = folder / "wrapper.toml"
    path.write_text("\n".join(lines) + "\n", encoding="utf8")
    return path


def run_radiopyo(args: argparse.Namespace, wrapper: Path, t_eval: np.ndarray | None):
    sys.path.insert(0, str(args.radiopyo))
    import radiopyo as rp

    logging.getLogger("radiopyo").setLevel(logging.ERROR)
    sim = rp.UnitCell.from_toml(wrapper)
    sim.change_LET(args.let)
    t_span = tuple(args.t_span)
    if args.mode == "run":
        return sim.run(t_span=t_span, t_eval=t_eval, max_step=args.max_step)
    if args.mode == "chunked":
        chunks = sim.prepare_chunked_run(
            t_span=t_span,
            max_step_size_on=args.max_step_on,
            max_step_size_off=args.max_step_off,
            max_step_size_final=args.max_step_final,
            points_per_chunk_on=args.points_on,
            points_per_chunk_off=args.points_off,
        )
    else:
        chunks = sim.prepare_single_pulse(
            t_span=t_span,
            max_step_size_on=args.max_step_on,
            max_step_size_trans=args.max_step_trans,
            max_step_size_final=args.max_step_final,
            points_per_chunk_on=args.points_on,
            points_per_chunk_off=args.points_off,
        )
    import collections

    return collections.deque(chunks, maxlen=1).pop()


def run_radkin(args: argparse.Namespace, wrapper: Path, out: Path, t_eval: np.ndarray | None) -> None:
    cmd = [str(args.radkin), "run", str(wrapper), "--t-span", repr(args.t_span[0]), repr(args.t_span[1]),
           "--mode", args.mode, "--let", repr(args.let), "--out", str(out),
           "--production", args.production, "--consumption", args.consumption]
    if t_eval is not None:
        f = out.parent / "t_eval.txt"
        f.write_text("\n".join(repr(float(t)) for t in t_eval), encoding="utf8")
        cmd += ["--t-eval", str(f)]
    for flag, value in (("--max-step", args.max_step), ("--max-step-on", args.max_step_on),
                        ("--max-step-off", args.max_step_off), ("--max-step-trans", args.max_step_trans),
                        ("--max-step-final", args.max_step_final)):
        if value is not None:
            cmd += [flag, repr(value)]
    if args.mode != "run":
        cmd += ["--points-on", str(args.points_on), "--points-off", str(args.points_off)]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        sys.exit(f"radkin failed ({done.returncode}):\n{done.stderr}")


def compare_frames(name: str, ref: pd.DataFrame, new: pd.DataFrame, rtol: float, atol: float) -> list[str]:
    errors = []
    if list(ref.columns) != list(new.columns):
        return [f"{name}: columns differ: {list(ref.columns)} vs {list(new.columns)}"]
    if len(ref) != len(new) or not np.allclose(ref.index, new.index, rtol=1e-12, atol=0):
        return [f"{name}: time grids differ ({len(ref)} vs {len(new)} points)"]
    worst = []
    for col in ref.columns:
        a, b = new[col].to_numpy(), ref[col].to_numpy()
        excess = np.abs(a - b) - (atol + rtol * np.abs(b))
        rel = np.max(np.abs(a - b) / np.maximum(np.abs(b), atol / rtol)) if len(b) else 0.0
        worst.append((rel, col))
        if np.any(excess > 0):
            i = int(np.argmax(excess))
            errors.append(f"{name} {col}: t={ref.index[i]:.6e} radkin={a[i]:.10e} radiopyo={b[i]:.10e}")
    worst.sort(reverse=True)
    print(f"  {name}: {len(ref.columns)} columns, {len(ref)} points, worst rel. diff "
          + ", ".join(f"{c}={r:.1e}" for r, c in worst[:3]))
    return errors


def eq_gross_floors(res, labels, rtol: float) -> dict[str, float]:
    """Floor of the net yield of each equilibrium pair: rtol x its gross
    (forward + backward) flux, since the net value is a difference of two
    large, nearly equal rates."""
    from scipy.integrate import trapezoid

    rates = res.reconstruct_reactions()
    out = {}
    for label in labels:
        if not res.env.reactions.is_eq_reaction(label):
            continue
        bi = res.env.reactions.get_eq_reaction_from(label)
        gross = rates[bi.forward.as_label()].abs() + rates[bi.backward.as_label()].abs()
        out[label] = rtol * 1e6 * float(trapezoid(gross.to_numpy(), x=gross.index.to_numpy()))
    return out


def compare_series(name: str, ref: pd.Series, new: pd.Series, rtol: float, floors=None) -> list[str]:
    ref, new = ref.astype(float), new.astype(float)
    if set(ref.index) != set(new.index):
        return [f"{name}: labels differ: {sorted(ref.index)} vs {sorted(new.index)}"]
    new = new[ref.index]
    floor = 1e-9 * float(np.max(np.abs(ref.to_numpy()))) if len(ref) else 0.0
    if floors:
        floor = pd.Series([max(floor, floors.get(k, 0.0)) for k in ref.index], index=ref.index)
    bad = np.abs(new - ref) > floor + rtol * np.abs(ref)
    rel = float(np.max(np.abs(new - ref) / np.maximum(np.abs(ref), floor / rtol + 1e-300))) if len(ref) else 0.0
    # rel: relative to max(|value|, floor/rtol), i.e. 1.0 x rtol at the pass limit
    print(f"  {name}: {len(ref)} values, worst rel. diff {rel:.1e}")
    return [f"{name} {k}: radkin={new[k]:.10e} radiopyo={ref[k]:.10e}" for k in ref.index[bad]]


def main() -> int:
    args = parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        wrapper = write_wrapper(args, tmp)
        t_eval = None
        if args.mode == "run":
            t0, t1 = args.t_span
            t_eval = np.geomspace(t0 + 1e-7 if t0 <= 0 else t0 + 1e-7, t1, args.t_eval_n)
        out = (args.keep or tmp / "radkin").resolve()
        run_radkin(args, wrapper, out, t_eval)
        res = run_radiopyo(args, wrapper, t_eval)

        errors: list[str] = []
        ref = res.to_pandas()
        new = pd.read_csv(out / "concentrations.csv", index_col=0)
        print(f"{args.config.name} mode={args.mode} beam={args.beam_type}:")
        errors += compare_frames("species", ref, new, args.rtol, 1e-10)
        ref_r = res.reconstruct_reactions()
        new_r = pd.read_csv(out / "reaction_rates.csv", index_col=0)
        peak = float(np.max(np.abs(ref_r.to_numpy()))) if ref_r.size else 0.0
        errors += compare_frames("rates", ref_r, new_r, args.rtol, 1e-7 * peak)
        new_i = pd.read_csv(out / "species_integrals.csv", index_col=0)["value"]
        errors += compare_series("integrate_species", res.integrate_species(), new_i, args.rtol)
        for sp in filter(None, args.production.split(",")):
            new_p = pd.read_csv(out / f"production_{sp}.csv", index_col=0)["value"]
            ref_p = res.integrate_production_of(sp)
            floors = eq_gross_floors(res, ref_p.index, args.rtol)
            errors += compare_series(f"production {sp}", ref_p, new_p, args.rtol, floors)
        for sp in filter(None, args.consumption.split(",")):
            new_c = pd.read_csv(out / f"consumption_{sp}.csv", index_col=0)["value"]
            ref_c = res.integrate_consumption_of(sp)
            floors = eq_gross_floors(res, ref_c.index, args.rtol)
            errors += compare_series(f"consumption {sp}", ref_c, new_c, args.rtol, floors)
        meta = json.loads((out / "run.json").read_text(encoding="utf8"))
        g_ref = {r.as_label(): r.ge() for r in res.env.reactions.radiolytic}
        for k, v in g_ref.items():
            if abs(meta["g_values"][k] - v) > 1e-12 * abs(v):
                errors.append(f"G-value {k}: radkin={meta['g_values'][k]} radiopyo={v}")
    if errors:
        print(f"FAIL ({len(errors)} difference(s)):")
        for e in errors[:25]:
            print("   ", e)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
