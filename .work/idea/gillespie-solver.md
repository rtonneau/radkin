---
status: grilled
scope: project
created: 2026-11-14
reviews: 0
promoted_to: null
---

## Problem

radkin solves only the deterministic mean-field rate equations (LSODA via rsoda). At low molecule counts (small voxels, early times, low dose), the mean-field description breaks down, and there is no independent solver to cross-check rsoda. Add a second solver based on the Gillespie algorithm that takes the same inputs as the LSODA path (SimEnv reactions and species, G-values, beams, `[voxel]` config).

## Decisions

- **Purpose: both.** It is validated against LSODA in the large-volume limit (ensemble mean → ODE solution), then used for small-volume stochastic physics.
- **Scope: unit cell and voxel grid.** Unit-cell SSA on the same SimEnv. On the voxel grid, an RDME: diffusion becomes hop events between neighbouring voxels (in place of Strang splitting), seeded by Geant4 molecules.
- **Algorithm:** exact SSA (direct / next-reaction method) as the reference, with tau-leaping as an opt-in speed-up.
- **Volume:** a new explicit config key gives the unit-cell volume used to convert mol/L to molecule counts. Voxel mode uses the voxel volume from `[voxel]`.
- **Beams:** each G-value × dose rate becomes a zero-order source reaction with a time-dependent propensity. Pulses are handled by SSA with time-varying propensities.
- **Output:** the ensemble of N trajectories (run in parallel with rayon, each with a reproducible per-trajectory seed) on the same time grid. The CSV layout matches LSODA's, plus std columns.
- **Validation:** at large volume, the ensemble mean lies within k·SEM of the LSODA result.
- **Selection:** a `solver = "lsoda" | "ssa"` config key, overridable by a CLI flag; the commands are otherwise unchanged.

## Open questions

- Which tau-leap variant to use (e.g. Cao–Gillespie–Petzold step selection), and how to handle the stiff fast equilibria (critical-reaction handling).
- How to convert rate constants to propensities for reactions between the same species (A + A: k vs 2k), matching radiopyo's rate convention.
- Which RNG crate and seeding scheme to use.
- How RDME interacts with voxel coarsening (ADR 0002 / `voxel-adaptive-coarsening`) and with the Geant4 seeds. Hop rates D/h² and the minimum voxel size at which RDME is still valid.

## Risks

- Fast equilibria (e.g. H+ + OH−) can make exact SSA spend almost all its events on them, making runs impractically slow.
- A wrong rate-convention factor would bias the comparison with LSODA without any visible error.
- RDME cost on large grids, and the loss of validity of reactions when voxels are too small.
- Statistical tests may be flaky; they need fixed seeds and sound tolerance bounds.

## Changelog

- 2026-11-14 — grill
