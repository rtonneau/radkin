---
bump: minor
floor: minor
---
- Voxel runs write a checkpoint (`checkpoint.npz` in the output directory) at `t_end`, and every `checkpoint_every` snapshot times when set in `[voxel.run]`.
- `radkin voxel CONFIG --resume CKPT --out DIR` resumes a voxel run from a checkpoint: crash recovery, extension to a later `t_end`, or new conditions (beam, solver, new seeds after the checkpoint time).
- Checkpoints are NumPy `.npz` archives readable with `np.load`.
