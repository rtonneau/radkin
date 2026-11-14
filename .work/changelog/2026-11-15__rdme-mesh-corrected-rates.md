---
bump: minor
floor: minor
---
- Stochastic voxel runs (`solver = "ssa"`) now correct bimolecular propensities for the voxel size by default (`[ssa] rate_correction = "hellander"`), so nm grids keep the macroscopic reaction rates; `rate_correction = "none"` gives the previous results. `snapshots.json` records which model ran.
- The RDME's critical-spacing warning now gives the right value for `2 A` reactions (it was half).
