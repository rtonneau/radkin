# Session: reduce-computational-burden

**Date:** 2026-11-14T13:04:58.948Z
**Status:** Grill phase complete

## Problem Statement

Reduce the computational cost of voxel mode (Phase 2, `radkin voxel`). The Phase-1 homogeneous solve was already optimized (~22%, session 2026-11-13, `docs/perf/README.md`); the per-voxel reaction work is where the cost multiplies.

## Context & Constraints

- **Current behavior:** `VoxelRun::advance_to` (`src/voxel/driver.rs`) uses step-doubling: one Strang step of `dt` plus two of `dt/2`, i.e. 3 Strang steps = 6 half-reactions per voxel per attempted step. Each half-reaction calls `solver::solve` cold per voxel (new Lsoda, Workspace, sparse pattern, order-1 start, `t_eval` output Vecs) via `reaction_step_background` (`src/voxel/reaction.rs`), parallel over rayon.
- **Pain point:** cold restarts and redundant solves dominate voxel runtime.
- **Dependencies:** keep radiopyo golden tests and Phase-1 results unchanged (the unit-cell path is not touched); existing voxel unit tests stay green; checkpoint/output formats unchanged; rsoda changes go to the rsoda repo and radkin bumps the pinned `rev`.
- **Tech stack:** Rust 2024, rsoda (LSODA), rayon (`parallel` feature).

## Success Metrics

- A voxel benchmark on the Geant4-seeded run (new `examples/voxel_bench.rs`, results recorded in `docs/perf/`) with baseline and per-lever numbers.
- Each lever is kept only if it cuts the plain median voxel run time by at least 5% on the benchmark (the docs/perf keep rule).
- Accuracy gate: against a reference run at tight control rtol (e.g. 1e-6), the new code at default tolerances has per-species domain totals and max per-voxel error no worse than the current code's error vs that reference.
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` clean.

## Architecture & Approach

Three levers, implemented and measured in this order:

1. **Lean per-voxel solve:** a dedicated integrate-to-t path for the reaction step that reuses per-thread Lsoda/workspace/sparse pattern, stores no t_eval output, and carries the last internal step size of each voxel as `first_step` of its next solve.
2. **Fuse half-reactions:** in the `dt/2` path, the two adjacent `dt/4` half-reactions (no diffusion between) become one `dt/2` solve: 6 to 5 voxel solves per attempted step.
3. **Embedded error estimate:** replace step-doubling with an estimate that reuses solves of the step (e.g. Strang vs a lower-order splitting on the same step), about 3 to about 1.5 Strang steps per accepted step. Recorded in an ADR.

## Assumptions & Trade-offs

- Voxel results may change within tolerance (not bit-identical); Phase-1/golden results must not change.
- Phase-1 solve optimizations are out of scope.
- Smarter background grouping (merging reaction solves of near-identical voxels) was not chosen.
- Adaptive grid coarsening (voxel merging as concentrations homogenise) is out of scope here: captured as idea `.work/idea/voxel-adaptive-coarsening.md` (two-way coarsen/refine, spread<tol merge with volume average, outputs/checkpoints prolonged to the fine grid) for a separate high-risk session.

## Open Questions

- Exact embedded estimator form (lever 3) is an implementation choice, constrained by the accuracy gate; documented in the ADR.
- Whether lever 1 needs an rsoda API change (re-init with a given first step / reuse of allocations).

## Notes

- Benchmark case: the Geant4-seeded setup used in session 2026-11-14__geant4-radkin-irt-syn-gvalues.
