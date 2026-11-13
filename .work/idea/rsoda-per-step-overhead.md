---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

Profile rsoda + radkin kernel per RHS evaluation and cut solver overhead

- Benchmark (tools/benchmark.py): radkin ~1.2 us per RHS evaluation incl. solver vs radiopyo ~37 us; 23-33x faster solve.
- Unknown split between kernel (derive/jacobian) and rsoda bookkeeping/LU; worth profiling before Phase 2 multiplies the cost per voxel.
- Ideas: avoid per-step allocations, reuse LU workspace, check Jacobian reuse heuristics vs ODEPACK.

## Changelog

- 2026-11-13 — capture — Profile rsoda + radkin kernel per RHS evaluation and cut solver overhead (radkin @ main)
