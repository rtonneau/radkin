# Session: rdme-reversible-mesh-correction

**Date:** 2026-11-15T11:50:50.802Z
**Status:** Grill phase complete

## Problem Statement

RDME mesh correction (`[ssa] rate_correction = "hellander"`) breaks detailed balance for reversible pairs, and the RDME voxel yields are not net of the background. In a seedless 2 um box at pH 7, about 480 H3O+ and OH- lose 8, 19 and 26 over 10 us at 20, 10 and 5 nm voxels, while the well-mixed unit-cell SSA keeps them constant. Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`, runs in `runs/_bg_rdme/`. Promoted from Idea `rdme-reversible-mesh-correction`.

## Context & Constraints

- **Current behavior:** each bimolecular propensity is scaled by 1/(1 - h*/h) on its own. The unimolecular direction of a `<->` pair is left unscaled. `yields.csv` subtracts only the t_start counts.
- **Pain point:** K_eq changes with voxel size, so the background drifts, and that drift ends up in the G-values.
- **Dependencies:** forward/backward pairs from the parser (`src/chemistry/reactions.rs` pair struct), single-grid and multi-grid RDME (ADR 0004, ADR 0005), the voxel SSA outputs. The unit-cell paths and the radiopyo golden tolerances must not change.
- **Tech stack:** Rust 2024, `src/ssa/` (propensity model from `SimMatrices`, rdme, rdme_multigrid), rayon.

## Success Metrics

- A seedless box at 20, 10 and 5 nm keeps H3O+/OH- within statistical noise of the unit-cell SSA over 10 us (`ssa_*` test, small size), on both the single grid and the multi-grid.
- The two directions of each declared `<->` pair share one factor. The pairs are logged, and unit tests on the radiopyo network cover them.
- Voxel SSA `yields.csv` is background-net through a seedless companion run, and its noise is reported. The single-track 10 us G-value report is rerun and compared with the deterministic voxel run.

## Architecture & Approach

- Pairs: only the forward/backward pairs declared with `<->`. Nothing is inferred from one-way reactions.
- Factor: the Hellander factor of the direction that is bimolecular once fixed species (H2O) are dropped, applied to both directions so K_eq is kept (Fange 2010; Hellander 2012; CRDME 2018 eq. 30). On the multi-grid, the factor comes from the level where the bimolecular direction runs. A pair that is bimolecular on both sides keeps its own factor per side, with a log line saying K_eq shifts.
- Background run: one seedless RDME run on the same grid and settings, with its own RNG seed derived from the run's seed. Its counts are subtracted from the seeded run. On by default for voxel SSA; `[ssa] background_run = false` turns it off. The output records the run and its noise.
- Docs: `docs/ssa.md`, ADR 0006, glossary terms "Reversible pair" and "Background run".

## Assumptions & Trade-offs

- The background run doubles runtime and memory. That cost is accepted, and the opt-out flag exists for it.
- One background run adds noise to small yields. Averaging several runs is not done here.
- Not done here: the sparse seed field memory fix (a separate Idea), and auto-detecting pairs among one-way reactions.

## Open Questions

- None blocking. Pairs that are bimolecular on both sides do not exist in the current data, and how to handle them may be revisited.

## Notes

Milestones from the Idea: (1) failing seedless `ssa_*` test; (2) pair detection and logging in `SimMatrices`; (3) shared factor on the single grid and the multi-grid; (4) background run and net yields; (5) G-value report rerun; (6) docs and ADR.
