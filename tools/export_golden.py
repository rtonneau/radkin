#!/usr/bin/env python3
"""
Export radiopyo's golden data (and a few extra reference values computed with
radiopyo) as JSON files read by the radkin Rust tests.

Run once, in the radiopyo conda env, with the radiopyo checkout to read:

    E:/miniconda3/envs/radiobio/python.exe -I tools/export_golden.py \
        --radiopyo D:/DEV/PYTHON/radiopyo --out tests/golden
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path

import numpy as np


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _npz(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def _rows(arrays: dict, prefix: str, kind: str) -> dict:
    return {
        f"{kind}_labels": [str(x) for x in arrays[f"{prefix}__{kind}_labels"]],
        kind: arrays[f"{prefix}__{kind}"].tolist(),
    }


def export_derive(g, out: Path) -> None:
    import radiopyo as rp  # noqa: F401
    from radiopyo.simulation.ode_functions_solver import jacobian

    golden = json.loads((g.DATA_DIR / "derive.json").read_text())
    for config in golden:
        cell = g.build_cell(config)
        labels = g._dyn_species_labels(cell)
        states = g._derive_states(cell)
        golden[config]["states"] = {k: {lb: float(v[lb]) for lb in labels} for k, v in states.items()}
        golden[config]["species_order"] = labels
        golden[config]["cst_species_order"] = [
            sp.label for sp in cell.env.species.sorted_iter_cst_species()
        ]
        reactions = cell.env.reactions
        golden[config]["reaction_str"] = (
            [str(r) for r in reactions.iter_all_kreactions()]
            + [str(r) for r in reactions.radiolytic]
            + [str(r) for r in reactions.enzymatic]
        )
        golden[config]["reaction_label"] = (
            [r.as_label() for r in reactions.iter_all_kreactions()]
            + [r.as_label() for r in reactions.radiolytic]
            + [r.as_label() for r in reactions.enzymatic]
        )
        golden[config]["k_values"] = [
            float(k) for _, _, k in reactions.iter_all_kreactions_and_kvalues()
        ]
        golden[config]["n_reactions"] = {
            "acid_base": len(reactions.ab_reaction),
            "equilibrium": len(reactions.eq_reaction),
            "k": len(reactions.k_reaction),
            "enzymatic": len(reactions.enzymatic),
            "radiolytic": len(reactions.radiolytic),
        }
        golden[config]["initial_values"] = [float(v) for v in cell.env.initial_values()]
        y = np.array([states["synthetic"][lb] for lb in labels])
        jac = jacobian(g.T_BEAM_ON, y, cell.beam, cell.mat)
        golden[config]["jacobian"] = {"synthetic__beam_on": jac.tolist()}
        y = np.array([states["synthetic"][lb] for lb in labels])
        y[labels.index("OH_r")] = -1e-3
        golden[config]["overshoot_state"] = {lb: float(v) for lb, v in zip(labels, y)}
        golden[config]["overshoot_dydt"] = dict(
            zip(labels, (float(v) for v in g._derive_adapter(cell, g.T_BEAM_ON, y.copy(), True)))
        )
        golden[config]["overshoot_jacobian"] = jacobian(g.T_BEAM_ON, y, cell.beam, cell.mat).tolist()
    (out / "derive.json").write_text(json.dumps(golden, indent=1))


def export_beams(g, out: Path) -> None:
    arrays = _npz(g.BEAMS_FILE)
    params = {
        "constant": {"type": "constant", "dose_rate": 1000.0, "max_dose": 5.0, "start_time": 0.2},
        "singlepulse": {"type": "singlepulse", "dose_rate": 1e6, "max_dose": 2.0, "start_time": 0.001},
        "pulsed": {"type": "pulsed", "dose_rate_peak": 1e6, "period": 0.01, "on_time": 2e-6,
                   "max_dose": 7.0, "start_time": 0.001},
        "pulsed_whole_pulses": {"type": "pulsed", "dose_rate_peak": 1e6, "period": 0.01,
                                "on_time": 2e-6, "max_dose": 6.0, "start_time": 0.001},
        "macropulsed": {"type": "macro_pulsed", "dose_rate_peak": 1e6, "period": 0.1,
                        "on_time": 0.001, "micro_period": 1e-4, "micro_on_time": 2e-6,
                        "max_dose": 50.0, "micro_start_time": 0.0, "start_time": 0.01},
    }
    out_d = {}
    for name, (_beam, windows) in g._beam_cases().items():
        out_d[name] = {
            "params": params[name],
            "scalars": arrays[f"{name}__scalars"].tolist(),
            "windows": [list(w) for w in windows],
            "matrix": [arrays[f"{name}__{i}__matrix"].tolist() for i in range(len(windows))],
            "intervals": [arrays[f"{name}__{i}__intervals"].tolist() for i in range(len(windows))],
        }
    (out / "beams.json").write_text(json.dumps(out_d, indent=1))


def export_runs(g, out: Path) -> None:
    from radiopyo.simulation.unit_cell import thinned_times

    arrays = _npz(g.RUNS_FILE)
    out_d = {}
    for prefix in ("labarbe_run", "labarbe_curve", "chunked_run"):
        entry = {
            "time": arrays[f"{prefix}__time"].tolist(),
            "n_points": arrays[f"{prefix}__n_points"].tolist(),
        }
        entry.update(_rows(arrays, prefix, "species"))
        entry.update(_rows(arrays, prefix, "rate"))
        out_d[prefix] = entry
    out_d["chunked_run"]["g_labels"] = [str(x) for x in arrays["chunked_run__g_labels"]]
    out_d["chunked_run"]["g_values"] = arrays["chunked_run__g_values"].tolist()
    out_d["labarbe_curve_t_eval"] = np.geomspace(g.RUN_T_SPAN[0], g.RUN_T_SPAN[1], g.N_CURVE).tolist()
    out_d["thinned_times"] = [
        {"args": [a, b, n], "times": thinned_times(a, b, n).tolist()}
        for a, b, n in ((0.0, 1.0, 5), (1e-7, 1e-3, 10), (2.0, 2.0, 4))
    ]
    # Short Labarbe run with integrals
    cell = g.build_cell("labarbe")
    t_eval = np.geomspace(1e-7, 1e-3, 50)
    res = cell.run(t_span=(1e-7, 1e-3), t_eval=t_eval)
    df = res.to_pandas()
    out_d["labarbe_short_integrals"] = {
        "t_eval": t_eval.tolist(),
        "time": df.index.to_list(),
        "species": {c: df[c].to_list() for c in df.columns},
        "integrate_species": {k: float(v) for k, v in res.integrate_species().items()},
        "production_H2O2": {k: float(v) for k, v in res.integrate_production_of("H2O2").items()},
        "consumption_O2": {k: float(v) for k, v in res.integrate_consumption_of("O2").items()},
        "burst_dose_rates": list(res.burst_dose_rates or []),
    }
    (out / "runs.json").write_text(json.dumps(out_d, indent=1))


def export_paper(radiopyo: Path, out: Path) -> None:
    pr = _load(radiopyo / "examples" / "tonneau_2025_pmb" / "paper_runs.py", "paper_runs")
    points = []
    for scenario, file in (("methods", "paper_tonneau_2025.npz"),
                           ("original", "paper_tonneau_2025_original.npz")):
        grid = _npz(radiopyo / "tests" / "golden" / "data" / file)
        sc = pr.SCENARIOS[scenario]
        for dr in (1.0, 1e3, 1e6):
            for po2 in (0.1, 2.0, 5.0, 21.0):
                i = pr.index_of(grid["dose_rates"], dr)
                j = pr.index_of(grid["po2"], po2)
                k = 2  # 30 Gy
                points.append({
                    "scenario": scenario,
                    "config": Path(sc.config).name,
                    "let": sc.let,
                    "dose_rate": dr,
                    "po2": po2,
                    "o2_initial_um": pr.o2_from_po2(po2),
                    "dose": float(grid["doses"][k]),
                    "t_end": pr.T_END,
                    "max_step_size_on": float(grid["doses"][k]) / dr / pr.ON_STEPS,
                    "max_step_size_final": pr.MAX_STEP_OFF,
                    "h2o2_final": float(grid["h2o2_final"][i, j, k]),
                    "o2_final": float(grid["o2_final"][i, j, k]),
                    "ph_min": float(grid["ph_min"][i, j, k]),
                })
    (out / "paper_points.json").write_text(json.dumps({"rtol": 1e-6, "points": points}, indent=1))


EDGE_CASES = [
    "H_r + H_r-> H2", "Fe++ + B -> C", "HO2-+H->X", "2O2 + 2OH_minus -> X",
    "A-> B - 2B + 3 C", "A -> -B + C", "A+B->2C+D", "-> A + B", "A + B -> C - D",
    "OO_r->x", "O2_r_minus -- SOD >> H2O2", "2H2O2 -- catalase >> O2 + 2 H2O",
    "A -- E >> B + C - D", "Ab2c -> d", "A + B -> ",
]


def export_parser(radiopyo: Path, out: Path) -> None:
    import tomllib

    from radiopyo.parser import RPARSER

    strings = []
    files = list((radiopyo / "radiopyo" / "data").glob("*.toml"))
    files += list((radiopyo / "examples" / "tonneau_2025_pmb" / "data").glob("*.toml"))
    files += list((radiopyo / "docs" / "issues" / "fenton").glob("*.toml"))
    for f in files:
        cfg = tomllib.loads(f.read_text(encoding="utf8"))
        for section in cfg.get("reactions", {}).values():
            if isinstance(section, list):
                strings += [r["reaction"] for r in section]
    strings += EDGE_CASES
    result = {}
    for s in dict.fromkeys(strings):
        gen = RPARSER.parse_reaction(s)
        result[s] = {
            "reactants": [[k, v[1], str(v[2])] for k, v in gen["reactants"].items()],
            "products": [[k, v[1], str(v[2])] for k, v in gen["products"].items()],
            "enzyme": None if gen["enzyme"] is None else gen["enzyme"].label,
        }
    (out / "parser_labels.json").write_text(json.dumps(result, indent=1))


def export_dense(out: Path) -> None:
    """scipy's LSODA dense output on real steps (Robertson), for src/dense.rs.

    Each record is the history as `LSODA._dense_output_impl` reads it
    (`rwork[20:]`, `iwork[13:15]`, `rwork[10:12]`) and the values of
    `solver.dense_output()` at points of the step just taken.
    """
    from scipy.integrate import LSODA

    def rober(_t, y):
        return [
            -0.04 * y[0] + 1.0e4 * y[1] * y[2],
            0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] ** 2,
            3.0e7 * y[1] ** 2,
        ]

    n = 3
    solver = LSODA(rober, 0.0, [1.0, 0.0, 0.0], 4.0e10, rtol=1e-4, atol=[1e-6, 1e-10, 1e-6])
    steps = []
    while solver.status == "running" and len(steps) < 400:
        solver.step()
        it = solver._lsoda_solver._integrator.iwork
        rw = solver._lsoda_solver._integrator.rwork
        nqu, nq = int(it[13]), int(it[14])
        ncol = max(nqu, nq) + 1
        dense = solver.dense_output()
        ts = [solver.t_old + f * (solver.t - solver.t_old) for f in (0.0, 0.25, 0.5, 0.9, 1.0)]
        steps.append({
            "yh": rw[20:20 + ncol * n].tolist(),
            "nyh": n,
            "n": n,
            "nq": nq,
            "nqu": nqu,
            "h": float(rw[11]),
            "hu": float(rw[10]),
            "tn": float(rw[12]),
            "t_step": float(solver.t),
            "t": ts,
            "y": [dense(t).tolist() for t in ts],
        })
    lowered = sum(s["nq"] < s["nqu"] for s in steps)
    assert lowered > 0, "no order decrease in the Robertson run"
    (out / "dense_output.json").write_text(json.dumps({"steps": steps}, indent=1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radiopyo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--only", choices=["dense"], help="export only this file")
    args = parser.parse_args()
    if args.only == "dense":
        args.out.mkdir(parents=True, exist_ok=True)
        export_dense(args.out)
        return
    sys.path.insert(0, str(args.radiopyo))
    import radiopyo  # noqa: F401

    logging.getLogger("radiopyo").setLevel(logging.ERROR)
    args.out.mkdir(parents=True, exist_ok=True)
    g = _load(args.radiopyo / "tests" / "golden" / "test_golden.py", "golden_scenarios")
    logging.getLogger("radiopyo").setLevel(logging.ERROR)
    export_parser(args.radiopyo, args.out)
    export_derive(g, args.out)
    export_beams(g, args.out)
    export_runs(g, args.out)
    export_paper(args.radiopyo, args.out)
    export_dense(args.out)
    for p in sorted(args.out.glob("*.json")):
        print(f"{p.name}: {p.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main()
