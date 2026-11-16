# Ticket 07: rdme-handover

**Status:** ✅ Done

## Local Test Result

`cargo test --release`: every test binary passes (33 `test result: ok`, no failure), including `irt_handover` (5 tests) and the existing `ssa_` tests unchanged. `cargo test --release irt_`: 48 passed.

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`.

## Review Notes

Base: `9965604cb87af8d5a53ac15f4778e700779bfd61`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- `radkin voxel` with an `[irt]` section and the ssa solver runs, per trajectory, the IRT stage to the handover time (`irt.t_end`), bins its molecules onto the fine grid and carries on with the RDME (`run_rdme_from_irt` in `src/ssa/voxel_sim.rs`). Trajectory `i` of the stage draws from stream `i` of `handover_seed(seed)`, apart from the RDME's own stream. Both pair rules work.
- Conservation, tested exactly for both rules: with `background_run = false` the net molecules of every record before the handover equal the particle counts of the same IRT stage run on its own, and at the handover the grid holds every one of them (`molecules_handed_over` equals the particle count times the trajectories).
- One time axis: `yields.csv` starts at the start of the IRT stage. Before the handover `N` is the particles plus the grid's Poisson draw of the initial solution (482 H3O+ and OH- at pH 7 in the test), so `N` is continuous and `G` is the particles per 100 eV.
- A seed at or after the handover time is refused ("outside the IRT stage").
- Without `[irt]` nothing changes: `run_rdme` takes the same path as before (`record` with no stage), `snapshots.json` gets no new key, and the existing `ssa_` tests pass unchanged.

Also checked: after a 100 ns handover the RDME (62.5 nm voxels) ends at 1 us within 25 % of the IRT stage carried on to 1 us, for e_aq, OH and H2O2.

Refused, with a message: `[irt]` with the deterministic solver; walls other than reflective; an `[irt.box]` that is not the grid; `trajectories` or `random_seed` in `[irt]` (they are those of `[ssa]`); a handover at or after `voxel.run.t_end`; a snapshot before the handover (the grid holds no field then); seeds left on the grid.

Choices:
- The grid starts at the handover time, from the Poisson draw of the initial solution; it does not simulate the background before the handover. The seedless background run (`[ssa] background_run`) does the same.
- `voxel.run.t_start` is not used in a run with an IRT stage: the run starts at the first seed, or at `irt.t_start`.
- The seeds are those of `[[irt.seed]]`, or of `[[voxel.seed]]` when `[irt]` has none; they are not binned onto the grid (`VoxelSimulation::from_toml_without_seeds`), which also avoids the dense seed field of the grid loader.
- The planned `src/irt/handover.rs` is not a file of its own: the handover lives next to the trajectory driver it extends (`Handover` in `voxel_sim.rs`), which keeps the RDME internals private.

Not in the ticket's file list: `src/ssa/mod.rs` (exports), `src/ssa/voxel_output.rs` (`irt_stage` in `snapshots.json`). `src/voxel/config.rs` did not need a change.

## Blockers / Challenges

None.

## Commits

- 180b4fe feat: IRT stage handed over to the RDME in radkin voxel (ticket 07)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
