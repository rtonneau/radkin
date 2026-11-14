---
status: promoted
scope: project
created: 2026-11-14
reviews: 0
promoted_to: radkin:2026-11-14__voxel-adaptive-coarsening
---

## Problem

Voxel mode spends reaction and diffusion work on every fine voxel even after concentrations have homogenised. Background-voxel dedup (`reaction_step_background`) already skips reaction solves for voxels equal to the uniform background, but smooth non-background regions still pay full reaction cost, and ADI diffusion always sweeps the full fine grid. Coarsening should cut reaction and diffusion cost equally.

- Came out of the grill of session 2026-11-14__reduce-computational-burden; run as its own high-risk session after that one's levers land.

## Decisions

- **Global uniform coarsening first.** The whole box changes level at once (2x2x2 -> 1 per level). Each level is a plain uniform grid, so ADI diffusion (exact only on uniform lines, where Ax/Ay/Az commute) and the background dedup keep working unchanged.
- **Local octree/AMR is a later, separate step**, only if runs show that localized features block global coarsening. It would need finite-volume face fluxes, coarse/fine flux matching and a sparse implicit solve, since ADI no longer applies.
- **Merge criterion:** the box coarsens only when every 2x2x2 block passes the per-species spread test (spread over the 8 children and their face neighbours below atol + rtol*|mean|, control tolerances times a factor). Merged value = volume average, which conserves mass exactly.
- **Alignment:** a `[voxel] max_coarsen_levels = L` option, with default 0 (off: current runs are unchanged). Config load rejects dims not divisible by 2^L with a clear error. Dims are never rounded silently.
- **Refine on re-injection only:** refine to level 0 before a new Geant4 seed or pulse is deposited (piecewise-constant prolongation), then let the grid re-coarsen. Otherwise it only coarsens, since diffusion only smooths.
- **I/O:** outputs and checkpoints stay on the fine grid (prolonged). On resume, the level is recomputed from the state, so the format doesn't change.
- **Acceptance:** species totals within rtol 1e-6 of a fine-grid run and fields within the merge tolerance, plus a measured speedup on a representative voxel case per the docs/perf keep rule. Otherwise it is not kept.

## Open questions

- Value of the tolerance factor in the merge criterion.
- Whether step-doubling / Strang dt control needs adjustment after a level change.
- Coarsening check cadence: every macro step or every N steps.

## Plan

1. Baseline: a representative voxel case that homogenises, with its fine-grid timing and results recorded in docs/perf.
2. Config: `[voxel] max_coarsen_levels` (default 0) and validation that dims divide by 2^L. Existing runs unchanged.
3. Restrict/prolong between levels plus the global spread test. Mass conserved to round-off (unit tests).
4. Driver runs diffusion and reactions on the current level, coarsening at macro-step boundaries. Settles the check cadence and the dt handling after a level change.
5. Refine to level 0 before re-injection; outputs and checkpoints prolonged to the fine grid; resume recomputes the level.
6. Acceptance: totals within rtol 1e-6 of the fine run, tolerance factor tuned, speedup measured against the keep rule. Decision recorded on keeping it and on whether the local-AMR follow-up is needed.

## Risks

- Localized features (one fresh track, a re-injection front) block global coarsening, so the gain may be small for many runs. That would trigger the local-AMR follow-up.
- Prolongation error at refine time (piecewise-constant data meets a fresh deposit).
- Smoothing bias if the tolerance factor is loose.

## Changelog

- 2026-11-14 — capture — from the reduce-computational-burden grill (radkin @ main)
- 2026-11-14 — grill
- 2026-11-14 — plan
- 2026-11-14 — promote
