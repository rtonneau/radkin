---
bump: minor
floor: minor
---
- Voxel SSA (RDME) runs with `rate_correction = "hellander"` keep the equilibrium of reversible pairs declared with `<->`: both directions now share the mesh-correction factor of their bimolecular direction, on the single grid and the multi-grid. A pH 7 background no longer drifts with voxel size.
- Voxel SSA yields are now net of the background: a seedless background run on the same grid is subtracted in the G-values (the `N_` columns stay gross). It is on by default; `[ssa] background_run = false` turns it off. It roughly doubles the runtime.
