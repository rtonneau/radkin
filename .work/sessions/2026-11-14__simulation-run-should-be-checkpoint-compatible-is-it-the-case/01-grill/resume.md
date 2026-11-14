# Session: Simulation run should be checkpoint compatible. Is it the case?

**Date:** 2026-11-14T11:25:43.563Z
**Status:** Grill phase complete

## Problem Statement

Voxel simulation runs (`radkin voxel`) cannot be checkpointed. A killed or crashed run loses all progress, a finished run cannot be extended to a later end time, and a run state cannot seed a new run with different conditions.

## Context & Constraints

- **Current behavior:** No checkpoint/resume code exists anywhere in `src/`. Every voxel run starts at t=0 from the config; `VoxelRun` (src/voxel/driver.rs) holds state, time and step control in memory only.
- **Pain point:** Long rayon voxel runs are not recoverable, extendable or restartable.
- **Dependencies:** Unit-cell (Phase-1) runs are untouched. Golden tests and radiopyo tolerances stay as they are. Species order follows radiopyo.
- **Tech stack:** Rust 2024, serde/serde_json (already present), new pure-Rust deps for NumPy `.npz` (ndarray-npy / zip). No C dependency (HDF5 rejected: libhdf5 too heavy).

## Success Metrics

- A voxel run interrupted at a checkpoint and resumed with `--resume` matches the uninterrupted run within run rtol (test).
- A finished run's final checkpoint can be extended to a later `t_end`; a checkpoint can be resumed with a changed beam/seeds/solver config (warnings logged).
- Resume with mismatching grid dims, spacing or species order fails with a `RadkinError`; a checkpoint file is never left half-written (temp + rename).
- `np.load` reads a checkpoint in Python; fmt, clippy and tests stay clean on Ubuntu and Windows.

## Architecture & Approach

- Scope: voxel runs only. Uses: crash recovery, extend a finished run, restart with new conditions.
- Format: NumPy `.npz` written in pure Rust: concentration arrays plus a JSON metadata entry (format version, t, grid dims/spacing, species order, current dt, IDs of injected seeds, config summary).
- Cadence: `[voxel] checkpoint_every = N` output times, plus always one at the end of the run. Atomic write (temp file, then rename).
- Resume: `radkin voxel <config> --resume <ckpt.npz> --out <new dir>`. Strict check on grid dims, spacing, species order (hard error); all other config differences allowed and logged as warnings.
- Seeds: checkpoint concentrations are the initial state; already-injected seeds are skipped; new seeds with times after t_ckpt are injected on top.
- Exactness: within run tolerance, not bit-identical. Step control restarts from the saved dt, with Rannacher on the first step.
- Outputs: a resumed run writes to a new output directory; the user stitches files.

## Assumptions & Trade-offs

- Not bit-identical: full step-control history is not saved.
- No unit-cell checkpointing. No HDF5/NetCDF. No automatic output stitching.

## Open Questions

None.

## Notes

The user asked whether runs are checkpoint compatible: they are not. HDF5 was first chosen, then replaced by `.npz` after discussing libhdf5's weight (a Rust rewrite was judged unrealistic).
