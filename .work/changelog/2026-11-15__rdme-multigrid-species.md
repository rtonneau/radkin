---
bump: minor
floor: minor
---
- Stochastic voxel runs (`solver = "ssa"`) accept `[ssa.species_levels]`, which puts a species on a coarser level of the voxel grid (e.g. `O2 = 2`: blocks of 4x4x4 voxels), where it hops and reacts. A numerous background species then costs far fewer hops: aerated O2 on level 2 of a 100 nm grid gave 14 times fewer hops with unchanged results, as long as the track does not deplete it locally. Without the table, results are unchanged.
- Reactions between species on different levels fire in the voxel of the coarsest reactant, and their products appear where the finest consumed molecule was. Global coarsening (`voxel.max_coarsen_levels`) and `tau_leap` work with species levels; `snapshots.json` records them.
- The run logs a suggested level for each background species (from its diffusion coefficient); it is never applied. `docs/ssa.md` explains the rules, the validation and when a coarse level biases the result.
