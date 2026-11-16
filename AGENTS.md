# AGENTS.md

radkin is a Rust port of radiopyo: homogeneous chemical kinetics of water radiolysis
(G-values, rate equations, beams), solved with the LSODA solver rsoda. Phase 1 reproduces
radiopyo's results. Phase 2 adds a voxelised 3D geometry with species diffusion, seeded by
Geant4-DNA molecules and coupled to the chemistry by Strang splitting (`src/voxel/`, ADR 0002).
The IRT stage (`src/irt/`, ADR 0007) runs the chemistry on the molecules as particles from
1 ps, with Geant4-DNA's reaction-time formulas, alone or before a handover to the voxel grid.
See `README.md` for usage and `docs/superpowers/specs/` for the Phase-1 design.

## Related repositories

Paths use `${DEV_DIR}`, the user-level dev root environment variable. Read its value
from the environment (`$env:DEV_DIR`, `$DEV_DIR`). If it is unset, stop and ask; never
guess a path.

| Repo | Local checkout | GitHub |
|---|---|---|
| rsoda (LSODA solver) | `${DEV_DIR}/RUST/rsoda` | https://github.com/rtonneau/rsoda |
| radiopyo (reference implementation) | `${DEV_DIR}/PYTHON/radiopyo` | https://github.com/rtonneau/radiopyo |

- radkin depends on rsoda as a git dependency pinned to a revision (`Cargo.toml`). A
  change needed in rsoda is made in the rsoda repo, pushed, and then radkin bumps the
  `rev`. Do not switch to a `path` dependency in committed code. A temporary
  `[patch]` for local testing is fine.
- radiopyo is the reference implementation. When radkin's behaviour is in doubt, read
  radiopyo's source and match it.

## Build and test

```bash
cargo fmt --check
cargo clippy --release --all-targets -- -D warnings
cargo test --release
```

CI (`.github/workflows/ci.yml`) runs these three checks on Ubuntu and Windows. Keep all of
them clean.

Comparison with radiopyo needs a Python environment where radiopyo is importable. Run
Python with `-I`.

- `python -I tools/compare_radiopyo.py ...`: compares one case.
- `pwsh tools/acceptance.ps1 -Python <python> -Radiopyo ${DEV_DIR}/PYTHON/radiopyo`: runs
  the 27-case acceptance matrix.
- `python -I tools/benchmark.py --repeat 5`: compares runtimes. The radkin side is
  `examples/bench.rs`.
- `python -I tools/breakdown.py --repeat 10`: the time breakdown of radkin's solves
  (kernel, rsoda linear solver, rest) on the same cases, no radiopyo needed. It uses the
  `timing` cargo feature. Results and the keep rule for optimizations are in `docs/perf/`.

Comparison of the IRT stage with Geant4-DNA needs no radiopyo, only a release build and the
Geant4 reference runs (IRT_syn, made with `${DEV_DIR}/GEANT4/SIM/radkin` and kept under
`.scratch/`; ten events are in `tests/data/irt/`).

- `python -I tools/compare_irt_geant4.py --g4 <run dir> --same-pass --events 2000`: G(t)
  of every species against Geant4-DNA IRT_syn, event by event, in free space, from the
  seeds the G(t) pass wrote itself. Results are in `docs/irt.md`.
- `python -I tools/compare_irt_geant4.py --events 100 --boundary outflow`: the same against
  the first reference, which ran in a 2 um water box that loses molecules at its faces
  (hence the outflow box) and took its seeds from a separate pass.
- `python -I tools/export_g4_chem.py`: writes `data/chem_g4.toml` from Geant4's reaction
  table dump.
- `cargo run --release --example irt_scaling -- --tracks 1,10,100`: run time against the
  number of overlapping tracks (`docs/perf/irt-scaling.md`).

## Layout

- `src/chemistry/`: reaction parser (a port of radiopyo's grammar), species, reactions,
  and `SimEnv`.
- `src/config/`: TOML loading, includes, and radiopyo's merge rules.
- `src/physics/`: timers, beams, and the beam collection.
- `src/gvalues.rs`: Boscolo LET G-values. This is an exact port of scipy's quadratic
  spline.
- `src/matrices.rs`, `src/kernel.rs`: rates, the right-hand side and the analytic
  Jacobian. They work on one cell at a time with buffers passed in by the caller, so they
  can be called per voxel in Phase 2.
- `src/solver.rs`: wraps rsoda to behave like scipy `solve_ivp(method="LSODA")`.
- `src/unit_cell.rs`, `src/result.rs`, `src/output.rs`, `src/bin/radkin.rs`: the run,
  chunked and single-pulse drivers, the results, and the CSV/JSON outputs and CLI.
- `src/voxel/`: Phase 2. Grid and state, `[voxel]` config, ADI diffusion, per-voxel reaction
  step (rayon), Strang driver with step-doubling, Geant4 seeds, simulation and outputs
  (`radkin voxel`).
- `src/ssa/`: stochastic solver (`solver = "ssa"`, ADR 0003, `docs/ssa.md`). Propensity
  model built from `SimMatrices`, exact SSA and tau-leaping on the unit cell, the ensemble
  and its outputs; on the voxel grid the RDME (hop events, coarsening, tau-leaping) and
  its multi-grid form with species levels (`rdme_multigrid*.rs`, ADR 0005). Its
  integration tests are named `ssa_*` so that `cargo test --release ssa` runs them all.
- `src/irt/`: IRT stage (`radkin irt`, `[irt]`, ADR 0007, `docs/irt.md`). Every molecule
  is a particle; reactions are independent reaction times per pair and per molecule.
  Species table and typed reaction model, sampling kernel (Geant4-DNA's formulas), the
  event loops of the two pair rules (`incremental.rs`, the default; `solver.rs`, full
  resampling),
  walls, seeds and outputs. The handover to the RDME is in `src/ssa/voxel_sim.rs`. Its
  integration tests are named `irt_*` so that `cargo test --release irt_` runs them all.
- `tests/golden/*.json`: radiopyo reference data exported by `tools/export_golden.py`.
- `tests/data/irt/`: Geant4-DNA reference of the IRT stage (reaction table dump, ten
  events).
- `data/`: radiopyo's config TOMLs and the Boscolo CSV; `chem_g4.toml`, the Geant4-DNA
  option3 chemistry written by `tools/export_g4_chem.py`.

## Conventions

- Results must stay within radiopyo's tolerances. The golden tests check the right-hand
  side, the Jacobian and the G-values at 1e-12, and runs at rtol 1e-6. Do not loosen a
  golden tolerance to make a test pass. Find the difference against radiopyo instead.
- Python float behaviour (`//`, `%`, `round`, numpy Horner evaluation, `geomspace`) is
  ported on purpose (`src/pyfloat.rs`, `src/poly.rs`). Keep these ports exact.
- The reaction index order and the species order follow radiopyo. Output CSV columns
  depend on them.
- Use Rust 2024 edition, `rustfmt.toml` (max width 100) and `thiserror` errors
  (`RadkinError`). The kernel does not allocate.
- Work goes on a branch with a PR to `main`.
- The idea backlog is in `.work/idea/`.
