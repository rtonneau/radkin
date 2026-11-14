# Implementation Plan

**Session:** Simulation run should be checkpoint compatible. Is it the case?
**Date:** 2026-11-14T11:33:08.417Z
**Estimated effort:** 1 day

## Strategy

Build the checkpoint file format first (pure data, testable alone), then wire writing and resuming into the voxel simulation, then expose it in the CLI with the end-to-end equivalence test and docs.

## Tickets Overview

- **Ticket 1:** `src/voxel/checkpoint.rs`: `Checkpoint` struct, atomic `.npz` write/read (arrays + JSON metadata), compatibility check against a grid/species.
- **Ticket 2:** `checkpoint_every` setting, periodic + final checkpoint writing, and `VoxelSimulation::resume` (initial state, dt, background, seed skipping, config-difference warnings).
- **Ticket 3:** `radkin voxel --resume <ckpt>` CLI, interrupted-vs-uninterrupted rtol test, extend and new-conditions tests, README docs.

## Sequencing Rationale

Ticket 2 needs the file format of ticket 1; ticket 3 needs the library API of ticket 2. Each leaves fmt/clippy/tests green.

## Risks & Mitigation

- **Risk:** npz crate versions conflict with existing deps or pull C code → **Mitigation:** use `zip` (deflate via pure-Rust miniz_oxide) and write the `.npy` headers by hand if ndarray-npy is awkward; check `cargo tree`.
- **Risk:** seeds double-injected or lost around t_ckpt → **Mitigation:** checkpoint records injected seeds (time + content fingerprint); seeds at time == t_ckpt are injected before the output callback, so they count as injected; tests cover it.
- **Risk:** resumed run drifts beyond rtol because step control restarts → **Mitigation:** save current dt and the background vector; set Rannacher pending on resume; the test compares at run rtol.

## Assumptions

- "Output times" are the `[voxel.run] snapshots` times; checkpoints are written from the `run_schedule` output callback.
- Totals/background totals before t_ckpt are not carried into the new output dir (the user stitches outputs).
- Windows and Ubuntu CI only need pure-Rust deps.
