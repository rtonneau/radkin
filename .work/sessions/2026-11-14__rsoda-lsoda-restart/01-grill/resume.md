# Session: rsoda-lsoda-restart

**Date:** 2026-11-14T17:16:33.185Z
**Status:** Grill phase complete

## Problem Statement

Each per-voxel reaction solve in radkin's voxel step builds a new rsoda `Lsoda`. That allocates a rsoda `State` and a dense n x n matrix the sparse path never uses. rsoda needs `Lsoda::restart(t0, y0, h0)` so one solver can be reset in place and reused across solves. An opt-in warm start that carries the method and order from the previous solve is in scope too.

## Context & Constraints

- **Current behavior:** `integrate_to` (`src/solver.rs`) builds a `ChemSystem<'a>` (borrowing `mat`, `source`, the scratch `Workspace` and `SparseJacobian`) and calls `Lsoda::new` for every voxel solve in `reaction_step_background` (`src/voxel/reaction.rs`). The first step comes from `CARRY_FACTOR` (0.01) times the carried step.
- **Pain point:** allocator contention. A local prototype of restart was bit-identical and gave about -7% at 16 threads on the voxel bench, with no gain at 1 thread (ticket 02 of session 2026-11-14__reduce-computational-burden).
- **Dependencies:** rsoda is a git dependency pinned to `rev`. The change is made and pushed in the rsoda repo first, then radkin bumps the `rev`. A `[patch]` is allowed for local testing only, never committed. Golden tolerances must not change.
- **Tech stack:** Rust 2024, rsoda (LSODA, dense and sparse linear solvers), rayon per-worker scratch.

## Success Metrics

- rsoda: `restart` with the default full reset gives the same results as `Lsoda::new` with the same options. A test compares them on several states and in both the dense and sparse Jacobian modes.
- radkin: fmt, clippy and test stay clean. The golden and acceptance results are unchanged, and the voxel results with warm start off stay within tolerance of the current ones.
- The restart change is kept only if a re-measured voxel bench (`docs/perf/` keep rule) shows a gain above noise. The warm start is kept only if it needs fewer RHS calls (nfev) and less wall time at equal accuracy against the reference run.

## Architecture & Approach

- **rsoda:** add `Lsoda::restart(t0, y0, h0)`. It keeps the owned system `S`, the `State` and the linear solver, and by default resets every counter, the history, the Nordsieck array, the method and the order the way `new()` does. An opt-in warm-start option keeps the method (Adams/BDF) and the order from the previous solve. The first step stays the caller's `h0`.
- **radkin:** bump the rsoda `rev`. Only the voxel reaction step uses `restart`: for each `reaction_step_background` call, every rayon worker holds one `ChemSystem` and one `Lsoda` (built in the per-worker `new_scratch` init) and restarts it for each voxel, the background solve included. `mat` is fixed within a call and the source is `None`, so the `ChemSystem<'a>` borrows need no `Arc` or ownership change. The chunked and single-pulse drivers keep `new()`.
- **Warm start in radkin:** a `[voxel]` config option, off by default. The first step still comes from `CARRY_FACTOR` times the carried step.

## Assumptions & Trade-offs

- The solver is rebuilt once per half-step per thread, not kept across steps. That is enough to remove the per-voxel allocations.
- We are not changing the non-voxel drivers: they make one solve per run, so there is nothing to gain.
- Warm start is off by default, so the golden and acceptance paths never use it.
- The bit-identical restart is a test target, while the external contract is "within tolerance".

## Open Questions

- rsoda has no history at a restart, so starting above order 1 may not be valid (the Nordsieck higher derivatives are unknown). If the prototype shows it fails or loses accuracy, the warm start shrinks to carrying only the method. That is decided during implementation, from measurements.

## Notes

- Promoted from Idea `rsoda-lsoda-restart` (`.work/idea/rsoda-lsoda-restart.md`).
- The gain is small and depends on the thread count, so it may vanish under bench noise or with another allocator.
- An incomplete reset (a field `new()` sets but `restart` forgets) could leak state silently between voxels. The `restart` vs `new()` test guards against this.
