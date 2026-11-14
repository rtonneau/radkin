---
status: captured
scope: project
created: 2026-11-14
reviews: 0
promoted_to: null
---

## Problem

Voxel step-doubling: take the first half-reaction of both the dt and the dt/2 paths from one solve (stop at t + dt/4, continue to t + dt/2), 4 reaction solves per attempt instead of 5

- Suggested in ticket 04 of session 2026-11-14__reduce-computational-burden after the embedded estimators were dropped (step-doubling's kept two-half-step solution is ~4x more accurate than one Strang step, so embedded pairs need about twice the steps).
- Same estimator; results no longer bit-identical (the dt path's first half-reaction passes through t + dt/4). Measure with examples/voxel_bench.rs and the accuracy gate.

## Changelog

- 2026-11-14 — capture — from ticket 04 (radkin @ perf/voxel-reaction-cost)
