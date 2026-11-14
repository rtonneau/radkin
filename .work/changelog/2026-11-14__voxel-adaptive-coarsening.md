---
bump: minor
floor: minor
---
- Voxel mode can coarsen the whole grid as concentrations homogenise: set `[voxel] max_coarsen_levels` (default 0, off) and optionally `coarsen_tol_factor` (default 0.1). Grid dims must be divisible by 2^max_coarsen_levels.
- Outputs and checkpoints stay on the fine grid; a resumed run recomputes its coarsening level from the checkpoint state.
- Resuming a voxel run now applies the Rannacher damping only when a seed was injected at the checkpoint time, so resumed runs match uninterrupted ones much more closely.
- `voxel_bench` gains `--max-coarsen-levels` and `--coarsen-tol-factor` and compares a coarsened run with the fine run.
