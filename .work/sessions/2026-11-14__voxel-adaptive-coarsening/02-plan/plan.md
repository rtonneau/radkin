# Implementation Plan

**Session:** voxel-adaptive-coarsening
**Date:** 2026-11-14T15:17:51.679Z
**Estimated effort:** 2-3 days

## Strategy

Measure first, then build bottom-up: a baseline homogenising case, the config surface (off by default), pure restrict/prolong/spread-test functions with conservation tests, driver integration (coarsen after accepted steps, refine on inject), fine-grid I/O and resume, then the acceptance run and the keep decision. Each ticket leaves `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings` and `cargo test --release` clean, and default runs (`max_coarsen_levels = 0`) unchanged.

## Tickets Overview

- **Ticket 1:** Baseline: a representative homogenising voxel case in `examples/voxel_bench.rs` with fine-grid timing recorded in `docs/perf/voxel.md`.
- **Ticket 2:** Config: `[voxel] max_coarsen_levels` (default 0) and `coarsen_tol_factor` (default 0.1), with dims validated against 2^L.
- **Ticket 3:** Coarsening primitives: restrict (volume average), prolong (piecewise constant) and the global spread test on uniform grids, with conservation unit tests.
- **Ticket 4:** Driver: the whole grid sits at one coarsening level; coarsen after each accepted splitting step, refine to level 0 in `inject`; dt proposal kept.
- **Ticket 5:** I/O: snapshots and checkpoints prolonged to the fine grid; resume recomputes the level; formats unchanged.
- **Ticket 6:** Acceptance: coarsened vs fine runs on the baseline case (totals within rtol 1e-6), tune `coarsen_tol_factor`, measure the speedup, record the keep decision and the local-AMR verdict in `docs/perf/voxel.md`.

## Sequencing Rationale

The baseline must exist before any change, so the speedup is measured against unchanged code. Config comes before the code that reads it. The primitives are pure and testable alone, the driver builds on them, I/O depends on the driver's level state, and acceptance needs everything.

## Risks & Mitigation

- **Risk:** Localized features block global coarsening, so the gain is small → **Mitigation:** the baseline case is chosen to homogenise; ticket 6 records how often coarsening is blocked and decides on the local-AMR follow-up.
- **Risk:** Mass drift across level changes → **Mitigation:** conservation unit tests at round-off (ticket 3) and domain-total checks in driver tests (ticket 4).
- **Risk:** Smoothing bias from a loose tolerance → **Mitigation:** default factor 0.1, tuned against the fine run in ticket 6.
- **Risk:** Prolongation error at refine time → **Mitigation:** refine only at injection, where the Rannacher start already damps the new deposit; compare against the fine run.

## Assumptions

- ADI diffusion and background dedup work unchanged on any uniform grid with doubled spacing (`VoxelGrid` with dims/2^L and spacing*2^L).
- Step-doubling sees only time error; the merge tolerance alone controls the spatial error of merging.
- Golden tests and radiopyo tolerances are untouched; nothing changes when coarsening is off.
