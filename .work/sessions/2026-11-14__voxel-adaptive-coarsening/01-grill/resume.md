# Session: voxel-adaptive-coarsening

**Date:** 2026-11-14T15:12:33.720Z
**Status:** Grill phase complete

## Problem Statement

Voxel mode spends reaction and diffusion work on every fine voxel even after concentrations have homogenised. Background-voxel dedup (`reaction_step_background`) skips reaction solves only for voxels equal to the uniform background; smooth non-background regions still pay the full reaction cost, and ADI diffusion always sweeps the full fine grid. Coarsening should cut reaction and diffusion cost equally.

## Context & Constraints

- **Current behavior:** Uniform fine voxel grid; Strang splitting step reaction(dt/2) -> ADI diffusion(dt) -> reaction(dt/2), dt adapted by step-doubling (time error only, max over voxels); `inject` adds seeds and sets the Rannacher start; outputs and checkpoints are on the fine grid.
- **Pain point:** Cost scales with the fine voxel count even when the field is smooth everywhere.
- **Dependencies:** ADI is exact only on a uniform grid (the 1D operators commute), so every level must stay uniform. Background dedup, checkpoint format and output formats stay unchanged. Results must stay within the run tolerances.
- **Tech stack:** Rust, rayon per-voxel reaction step, rsoda, `src/voxel/` (grid, diffusion, driver, sim, checkpoint, output, config).

## Success Metrics

- With `max_coarsen_levels = 0` (default), existing runs and tests are unchanged.
- Restriction/prolongation conserves the domain totals of every species to round-off (unit tests).
- On a representative homogenising voxel case, species totals stay within rtol 1e-6 of the fine-grid run and fields within the merge tolerance.
- That case gets measurably faster, per the docs/perf keep rule; otherwise the feature is not kept.

## Architecture & Approach

Global uniform coarsening. The whole grid sits at one **coarsening level** L (each coarse voxel = 2^L x 2^L x 2^L fine voxels). After each accepted splitting step, if `L < max_coarsen_levels` and every 2x2x2 block passes the per-species spread test (spread over the 8 children and their face neighbours below `coarsen_tol_factor * (atol + rtol*|mean|)`), the grid drops one level: merged value = volume average (mass conserved exactly), the grid spacing doubles, ADI and background dedup run unchanged on the coarser uniform grid. The dt proposal is kept after a level change (merging only smooths, and implicit ADI has no CFL limit).

Refine: `inject` prolongs piecewise-constant to level 0 before each deposit (seed or pulse), then the grid may re-coarsen. Snapshots and checkpoints are prolonged to the fine grid; resume recomputes the level from the state, so no format change.

Config: `[voxel] max_coarsen_levels` (default 0 = off) and `[voxel] coarsen_tol_factor` (default 0.1). Config load rejects dims not divisible by 2^L with a clear error; dims are never rounded silently.

Milestones: (1) baseline homogenising case with fine-grid timing in docs/perf; (2) config + validation; (3) restrict/prolong + global spread test with conservation tests; (4) driver coarsens after accepted steps on the current level; (5) refine on injection, fine-grid outputs/checkpoints, resume; (6) acceptance and keep decision, including whether the local-AMR follow-up is needed.

## Assumptions & Trade-offs

- Assumes many runs homogenise everywhere at a similar rate; one localized feature (fresh track, re-injection front) blocks global coarsening, so the gain may be small for those runs.
- Not doing now: local octree/AMR (needs FV face fluxes, coarse/fine flux matching and a sparse implicit solve). A separate follow-up only if the acceptance milestone shows localized features block coarsening.
- Not doing: per-axis (anisotropic) levels, auto-rounding dims, refine indicators between injections (diffusion only smooths).
- Step-doubling does not see the spatial error of merging; only the merge tolerance controls it.

## Open Questions

None blocking. The value of `coarsen_tol_factor` (default 0.1) is tuned in the acceptance milestone.

## Notes

- Came from the Idea `voxel-adaptive-coarsening` (grilled, then planned 2026-11-14), originally from the reduce-computational-burden grill.
- Risks: prolongation error at refine time (piecewise-constant data meets a fresh deposit); smoothing bias if the tolerance factor is loose.
- Glossary: added **Coarsening level**.
