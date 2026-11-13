#!/usr/bin/env python3
"""
Time breakdown of radkin's solves on the benchmark matrix (no radiopyo needed).

For every case of tools/benchmark.py, dense and sparse, this script runs
examples/bench.rs twice:
  - plain build: the median solve time, with no timer in the code;
  - `timing` build (target/timing): `--breakdown`, the time spent in kernel
    RHS calls, kernel Jacobian calls, rsoda's linear-solver `prepare`
    (Jacobian evaluation + assembly + factorization) and `solve` calls.
Shares are given relative to the instrumented median solve time:
  rhs      kernel derive
  jac      kernel Jacobian (dense or sparse)
  factor   prepare - jac (assembly and LU factorization of P = I - h*el0*J)
  lin      linear solves with the factors
  rest     everything else: rsoda step bookkeeping and radkin's output storage
"timer" is the instrumented median over the plain median (timer overhead).

    python -I tools/breakdown.py --repeat 10 --out .scratch/breakdown.json
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from benchmark import CASES, t_eval_for, wrapper_for  # noqa: E402

EXE = "bench.exe" if sys.platform == "win32" else "bench"


def build(target_dir: Path | None, features: list[str]) -> Path:
    cmd = ["cargo", "build", "--release", "--quiet", "--example", "bench", *features]
    if target_dir is not None:
        cmd += ["--target-dir", str(target_dir)]
    if subprocess.run(cmd, cwd=ROOT).returncode != 0:
        sys.exit("cargo build failed")
    return (target_dir or ROOT / "target") / "release" / "examples" / EXE


def run_bench(exe: Path, case: dict, wrapper: Path, t_eval, repeat: int, folder: Path,
              extra: list[str]) -> dict:
    cmd = [str(exe), str(wrapper), "--t-span", *map(repr, case["t_span"]),
           "--mode", case["mode"], "--repeat", str(repeat), "--jacobian", "both", *extra]
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


def summarize(plain: dict, timed: dict) -> dict:
    med = statistics.median
    b = timed["breakdown"]
    comp = {k: med(r[k] for r in b) for k in ("t_rhs", "t_jac", "t_prepare", "t_solve")}
    total = med(timed["solve_s"])
    parts = dict(rhs=comp["t_rhs"], jac=comp["t_jac"], factor=comp["t_prepare"] - comp["t_jac"],
                 lin=comp["t_solve"])
    parts["rest"] = total - sum(parts.values())
    last = b[-1]
    return dict(
        solve_s=med(plain["solve_s"]), timed_s=total, nfev=plain["nfev"], nst=plain["nst"],
        njev=last["n_jac"], n_solve=last["n_solve"],
        per_rhs_us=med(plain["solve_s"]) / plain["nfev"] * 1e6,
        share={k: v / total for k, v in parts.items()}, seconds=parts,
        timer=total / med(plain["solve_s"]))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--repeat", type=int, default=10)
    p.add_argument("--only", help="run only cases whose name contains this text")
    p.add_argument("--out", type=Path, help="write the rows to this JSON file")
    args = p.parse_args()

    plain_exe = build(None, [])
    timed_exe = build(ROOT / "target" / "timing", ["--features", "timing"])
    rows = []
    for case in CASES:
        if args.only and args.only not in case["name"]:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            wrapper = wrapper_for(case, tmp)
            t_eval = t_eval_for(case)
            plain = run_bench(plain_exe, case, wrapper, t_eval, args.repeat, tmp, [])
            timed = run_bench(timed_exe, case, wrapper, t_eval, args.repeat, tmp, ["--breakdown"])
        for jac, pl, ti in (("dense", plain, timed), ("sparse", plain["sparse"], timed["sparse"])):
            row = dict(case=case["name"], jacobian=jac, **summarize(pl, ti))
            rows.append(row)
            print(f"{case['name']:<38} {jac:<6} solve {row['solve_s'] * 1e3:8.3f} ms", flush=True)

    print()
    print("| case | jac | solve [ms] | nfev | nst | njev | per-RHS [us] "
          "| rhs | jac | factor | lin | rest | timer |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        s = r["share"]
        print(f"| {r['case']} | {r['jacobian']} | {r['solve_s'] * 1e3:.3f} | {r['nfev']} | "
              f"{r['nst']} | {r['njev']} | {r['per_rhs_us']:.2f} | "
              + " | ".join(f"{100 * s[k]:.0f}%" for k in ("rhs", "jac", "factor", "lin", "rest"))
              + f" | x{r['timer']:.2f} |")
    if args.out:
        args.out.write_text(json.dumps(dict(repeat=args.repeat, rows=rows), indent=1),
                            encoding="utf8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
