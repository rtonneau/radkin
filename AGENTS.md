# AGENTS.md

radkin is a Rust port of radiopyo: homogeneous chemical kinetics of water radiolysis
(G-values, rate equations, beams), solved with the LSODA solver rsoda. Phase 1 reproduces
radiopyo's results. Phase 2 adds a voxelised 3D geometry with species diffusion, seeded by
Geant4-DNA molecules and coupled to the chemistry by Strang splitting (`src/voxel/`, ADR 0002).
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
- `tests/golden/*.json`: radiopyo reference data exported by `tools/export_golden.py`.
- `data/`: radiopyo's config TOMLs and the Boscolo CSV.

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
