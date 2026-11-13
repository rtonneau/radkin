#!/usr/bin/env python3
"""
Runtime comparison of radkin and radiopyo on the same cases.

Both codes are timed in-process, on the same wrapper TOML as
compare_radiopyo.py (config + beam), for --repeat repeats:
  - load:  build the simulation from the TOML (parsing, includes, matrices,
           G-values) — radiopyo `UnitCell.from_toml` (+ change_LET),
           radkin `UnitCell::from_toml` (+ change_let);
  - solve: the integration — `run`, `prepare_chunked_run(...)` or
           `prepare_single_pulse(...)` fully consumed / radkin's equivalents.
No output file is written by either code. radiopyo gets one untimed warm-up
run (lazy imports, caches); medians are reported.

radkin side: `cargo run --release --example bench` (examples/bench.rs).

    python -I tools/benchmark.py --radiopyo D:/DEV/PYTHON/radiopyo --repeat 5 \
        --out bench.json
"""
from __future__ import annotations

import argparse
import collections
import json
import logging
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from compare_radiopyo import write_wrapper  # noqa: E402

# name, config, mode, t_span, beam (type, params) or None, options, extra
CASES = [
    dict(name="labarbe run 1e-7..10 s", config="data/config_Labarbe.toml", mode="run",
         t_span=(1e-7, 10.0), beam=None, opts={}),
    dict(name="labarbe run, 200 t_eval", config="data/config_Labarbe.toml", mode="run",
         t_span=(1e-7, 10.0), beam=None, opts={}, t_eval_n=200),
    dict(name="water constant run 0..1 s", config="data/advanced_water_radiolysis.toml", mode="run",
         t_span=(0.0, 1.0), beam=("constant", {"dose_rate": 1e3, "max_dose": 10}), opts={}, let=0.201),
    dict(name="water single pulse chunked", config="data/advanced_water_radiolysis.toml",
         mode="chunked", t_span=(0.0, 1.0), beam=("singlepulse", {"dose_rate": 1e6, "max_dose": 10}),
         opts={"max_step_on": 1e-7, "max_step_off": 0.05}, let=0.201),
    dict(name="water pulsed 100 Hz chunked 0..1 s", config="data/advanced_water_radiolysis.toml",
         mode="chunked", t_span=(0.0, 1.0),
         beam=("pulsed", {"dose_rate_peak": 1e6, "period": 0.01, "on_time": 2e-6, "max_dose": 100}),
         opts={"max_step_on": 1e-7, "max_step_off": 1e-3}, let=0.201),
    dict(name="water macro-pulsed chunked 0..0.1 s", config="data/advanced_water_radiolysis.toml",
         mode="chunked", t_span=(0.0, 0.1),
         beam=("macro_pulsed", {"dose_rate_peak": 1e6, "period": 0.05, "on_time": 0.001,
                                "micro_period": 1e-4, "micro_on_time": 2e-6, "max_dose": 10}),
         opts={"max_step_on": 1e-7, "max_step_off": 1e-4}, let=0.201),
    dict(name="labarbe pulsed run 0..0.015 s", config="data/config_Labarbe.toml", mode="run",
         t_span=(0.0, 0.015),
         beam=("pulsed", {"dose_rate_peak": 1e6, "period": 0.01, "on_time": 2e-6, "max_dose": 10}),
         opts={"max_step": 1e-6}, let=0.201),
    dict(name="paper single pulse 0..1 s", config="data/CONFIG_PAPER_01.toml", mode="single-pulse",
         t_span=(0.0, 1.0), beam=("singlepulse", {"dose_rate": 1e6, "max_dose": 10}),
         opts={"max_step_on": 1e-7, "max_step_trans": 1e-5, "max_step_final": 0.05},
         let=0.56, variable_gvalues=True),
    dict(name="paper pulsed chunked 0..0.1 s", config="data/CONFIG_PAPER_01.toml", mode="chunked",
         t_span=(0.0, 0.1),
         beam=("pulsed", {"dose_rate_peak": 1e6, "period": 0.01, "on_time": 2e-6, "max_dose": 10}),
         opts={"max_step_on": 1e-7, "max_step_off": 1e-3}, let=0.56, variable_gvalues=True),
]


def wrapper_for(case: dict, folder: Path) -> Path:
    beam_type, params = case["beam"] if case["beam"] else (None, {})
    ns = argparse.Namespace(
        config=ROOT / case["config"], variable_gvalues=case.get("variable_gvalues", False),
        initial=[], beam_type=beam_type, beam_param=[f"{k}={v}" for k, v in params.items()],
        let=case.get("let", 0.201))
    return write_wrapper(ns, folder)


def t_eval_for(case: dict) -> np.ndarray | None:
    n = case.get("t_eval_n")
    if not n:
        return None
    t0, t1 = case["t_span"]
    return np.geomspace(max(t0, 0.0) + 1e-7, t1, n)


def radiopyo_once(rp, case: dict, wrapper: Path, t_eval):
    t0 = time.perf_counter()
    sim = rp.UnitCell.from_toml(wrapper)
    if "let" in case:
        sim.change_LET(case["let"])
    t1 = time.perf_counter()
    o = case["opts"]
    if case["mode"] == "run":
        res = sim.run(t_span=case["t_span"], t_eval=t_eval, max_step=o.get("max_step"))
    else:
        kw = dict(t_span=case["t_span"], max_step_size_on=o.get("max_step_on"),
                  max_step_size_final=o.get("max_step_final"))
        if case["mode"] == "chunked":
            chunks = sim.prepare_chunked_run(max_step_size_off=o.get("max_step_off"), **kw)
        else:
            chunks = sim.prepare_single_pulse(max_step_size_trans=o.get("max_step_trans"), **kw)
        res = collections.deque(chunks, maxlen=1).pop()
    t2 = time.perf_counter()
    return t1 - t0, t2 - t1, res


def bench_radiopyo(rp, case, wrapper, t_eval, repeat):
    radiopyo_once(rp, case, wrapper, t_eval)  # warm-up
    load, solve = [], []
    for _ in range(repeat):
        lo, so, res = radiopyo_once(rp, case, wrapper, t_eval)
        load.append(lo)
        solve.append(so)
    return dict(load_s=load, solve_s=solve, nfev=int(sum(res.nfev)),
                n_points=int(len(res.t)), n_chunks=len(res.n_points))


def bench_radkin(exe: Path, case, wrapper, t_eval, repeat, folder: Path):
    cmd = [str(exe), str(wrapper), "--t-span", *map(repr, case["t_span"]),
           "--mode", case["mode"], "--repeat", str(repeat)]
    if "let" in case:
        cmd += ["--let", repr(case["let"])]
    for k, v in case["opts"].items():
        cmd += ["--" + k.replace("_", "-"), repr(v)]
    if t_eval is not None:
        f = folder / "t_eval.txt"
        f.write_text("\n".join(repr(float(t)) for t in t_eval), encoding="utf8")
        cmd += ["--t-eval", str(f)]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        sys.exit(f"radkin bench failed on {case['name']}:\n{done.stderr}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--radiopyo", type=Path, default=Path("D:/DEV/PYTHON/radiopyo"))
    p.add_argument("--repeat", type=int, default=5)
    p.add_argument("--only", help="run only cases whose name contains this text")
    p.add_argument("--out", type=Path, help="write all timings to this JSON file")
    args = p.parse_args()

    build = subprocess.run(["cargo", "build", "--release", "--quiet", "--example", "bench"], cwd=ROOT)
    if build.returncode != 0:
        return 1
    exe = ROOT / "target" / "release" / "examples" / ("bench.exe" if sys.platform == "win32" else "bench")

    sys.path.insert(0, str(args.radiopyo))
    import radiopyo as rp

    logging.getLogger("radiopyo").setLevel(logging.ERROR)
    med = statistics.median
    rows = []
    for case in CASES:
        if args.only and args.only not in case["name"]:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            wrapper = wrapper_for(case, tmp)
            t_eval = t_eval_for(case)
            rk = bench_radkin(exe, case, wrapper, t_eval, args.repeat, tmp)
            py = bench_radiopyo(rp, case, wrapper, t_eval, args.repeat)
        row = dict(case=case["name"], mode=case["mode"], radkin=rk, radiopyo=py,
                   speedup_solve=med(py["solve_s"]) / med(rk["solve_s"]),
                   speedup_total=(med(py["load_s"]) + med(py["solve_s"]))
                   / (med(rk["load_s"]) + med(rk["solve_s"])))
        rows.append(row)
        print(f"{case['name']:<38} solve radiopyo {med(py['solve_s']):9.4f} s  radkin "
              f"{med(rk['solve_s']):9.5f} s  x{row['speedup_solve']:6.1f}   "
              f"nfev {py['nfev']}/{rk['nfev']}", flush=True)

    print()
    print("| case | radiopyo load [s] | radkin load [s] | radiopyo solve [s] | radkin solve [s] "
          "| solve speed-up | total speed-up | nfev radiopyo / radkin |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        py, rk = r["radiopyo"], r["radkin"]
        print(f"| {r['case']} | {med(py['load_s']):.4f} | {med(rk['load_s']):.5f} | "
              f"{med(py['solve_s']):.4f} | {med(rk['solve_s']):.5f} | "
              f"{r['speedup_solve']:.0f}x | {r['speedup_total']:.0f}x | {py['nfev']} / {rk['nfev']} |")
    if args.out:
        args.out.write_text(json.dumps(dict(repeat=args.repeat, rows=rows), indent=1), encoding="utf8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
