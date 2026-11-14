---
status: promoted
scope: project
created: 2026-11-14
reviews: 0
promoted_to: radkin:2026-11-14__rsoda-lsoda-restart
---

## Problem

Each per-voxel reaction solve in radkin's voxel step builds a new rsoda `Lsoda`. That allocates a rsoda `State` and a dense n x n matrix the sparse path never uses. rsoda needs `Lsoda::restart(t0, y0, h0)` so one solver can be reset in place and reused by every solve on a thread.

- Found in ticket 02 of session 2026-11-14__reduce-computational-burden. A local prototype gave bit-identical results and about -7% at 16 threads on the voxel bench, with no gain at 1 thread (allocator contention).
- The work lands in rsoda. radkin then bumps the pinned `rev` and uses `restart` from the per-thread `SolveScratch`.

## Decisions

- **API: `restart(t0, y0, h0)` keeps the owned system `S`.** rsoda keeps its `State` and linear solver across restarts. On the radkin side, `ChemSystem` becomes long-lived per thread: it lives in or next to `SolveScratch` and owns its `Workspace`. The per-solve inputs (source, matrices) are swapped through `system_mut` or a setter. Rejected: `restart_with(system, ...)` and a detachable `State` with `from_parts`/`into_parts`.
- **Contract: within tolerance, not bit-identical.** By default, `restart` does a full reset equivalent to `new()` with the same options: counters, history, Nordsieck array, method and order. An **opt-in warm-start flag** keeps the method (Adams/BDF) and order hints from the previous solve. Golden paths use the default.
- **Keep rule:** re-measure on the voxel bench after the rev bump (`docs/perf/` keep rule). Drop the radkin change if the gain is within noise.
- rsoda stays a git dependency pinned to a revision. A `[patch]` is used only for local testing.

## Open questions

- What exactly does the warm-start flag keep (method, order, step only)? Does it beat or replace `CARRY_FACTOR` (0.01) for the first step?
- Can `ChemSystem`'s borrows (`mat`, `source`, `ws`) be restructured so one system lives across solves without fighting the `Lsoda<'a, S>` lifetime? Or does it need to own or `Arc` its matrices?
- Should the warm-start flag be measured in this Idea, or split out as a follow-up?

## Risks

- The gain is small (~7%) and shows only at high thread counts, so it may vanish under bench noise or with another allocator.
- An incomplete reset (a field `new()` sets but `restart` forgets) could leak state silently between voxels. It needs a test comparing `restart` to `new()` on several states and both Jacobian modes.
- Making `ChemSystem` long-lived may ripple through `solver.rs` and the chunked/single-pulse drivers.
- The warm-start option could break the rtol 1e-6 golden runs if it were ever enabled on those paths.

## Changelog

- 2026-11-14 — capture — from ticket 02 review (radkin @ perf/voxel-reaction-cost)
- 2026-11-14 — grill
- 2026-11-14 — promote
