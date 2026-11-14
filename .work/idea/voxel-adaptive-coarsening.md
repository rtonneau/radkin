---
status: captured
scope: project
created: 2026-11-14
reviews: 0
promoted_to: null
---

## Problem

Voxel mode: adaptive grid coarsening (octree/AMR) so the cell count drops as concentrations homogenise

- Decided in the grill of session 2026-11-14__reduce-computational-burden; to be run as its own high-risk session after that one's levers land.
- Two-way: coarsen 2x2x2 blocks and refine again (re-injection, sharp fronts).
- Merge criterion: per-species spread over the 8 children and their face neighbours below atol + rtol*|mean| (control tolerances times a factor); merged value = volume average, which conserves mass exactly.
- Outputs and checkpoints stay on the fine grid (merged cells prolonged piecewise-constant); resume rebuilds the tree from the state. No format change.
- Open: non-uniform diffusion (ADI no longer applies directly), flux matching at coarse/fine faces, refine indicator.

## Changelog

- 2026-11-14 — capture — from the reduce-computational-burden grill (radkin @ main)
