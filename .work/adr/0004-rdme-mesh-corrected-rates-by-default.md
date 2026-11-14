# Stochastic voxel runs correct bimolecular propensities for the voxel size by default

_ADR 0004, recorded 2026-11-15 (session 2026-11-15__rdme-mesh-corrected-rates)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

With the macroscopic constant the RDME underestimates bimolecular encounters once the voxel spacing approaches the critical spacing `h* = 0.2527 k / D` (a factor of 1.2 to 3 on radical recombinations at 6 nm). `[ssa] rate_correction` therefore defaults to `"hellander"`: each bimolecular propensity on a cubic grid level of spacing `h > h*` is scaled by `1 / (1 - h*/h)` (Hellander & Hellander 2020, eq. 6, with `1/k = 1/k_r + 1/(4 pi sigma D)`). It needs no reaction radius. Reactions with `h <= h*` stay uncorrected and get a warning, and anisotropic grids fall back to uncorrected with a warning unless the key is set explicitly. Every stochastic voxel result changes; `rate_correction = "none"` reproduces the earlier model.

Considered: opt-in (no result changes, but the biased model stays the default); refusing runs with `h <= h*` (rules out the bundled chemistry at nm voxels, where `e_aq + O_minus + H2O` has `h* = 75 nm`).
