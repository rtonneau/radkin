---
status: captured
scope: project
created: 2026-11-14
reviews: 0
promoted_to: null
---

## Problem

rsoda: `Lsoda::restart(t0, y0, h0)` to reset a solver in place, so radkin's voxel reaction step stops allocating rsoda `State` and the unused dense n x n matrix per voxel solve

- Found in ticket 02 of session 2026-11-14__reduce-computational-burden; a local prototype was bit-identical and gave about -7% at 16 threads on the voxel bench, no gain at 1 thread.
- Work lands in rsoda, then radkin bumps the pinned rev and uses it from the per-thread `SolveScratch`.

## Changelog

- 2026-11-14 — capture — from ticket 02 review (radkin @ perf/voxel-reaction-cost)
