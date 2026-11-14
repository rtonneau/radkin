# Ticket 05: io-resume

**Status:** ✅ Done

## Local Test Result

Base commit 5a068121276adaffc60732f406b1ae05c479b88d (clean tree).

- `cargo test --release voxel::`: 58 passed (lib), 0 failed. Output in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/evidence/05-voxel-tests.txt`.
- `cargo fmt --check`: clean.
- `cargo clippy --release --all-targets -- -D warnings`: clean.
- `cargo test --release` (full): all suites pass (77 lib tests, `voxel_checkpoint` 6 passed, others unchanged). Output in `evidence/05-full-tests.txt`.

## Review Notes

Base commit: 5a068121276adaffc60732f406b1ae05c479b88d. No review rounds yet.

- AC1 (fine-grid snapshots/CSV/JSON): ticket 04 already routes snapshots through `VoxelRun::fine_state()`. New test `snapshots_outputs_and_checkpoints_stay_on_the_fine_grid_when_coarsened` (tests/voxel_checkpoint.rs) runs the example at level 2 and checks: the last snapshot has fine n_voxels and is piecewise constant over 4x4x4 blocks; `write_voxel_dir` gives the same `snapshots.json` keys/dims/file list, same `.npy` file names, sizes and headers as an uncoarsened run; the last `yields.csv` row equals the fine snapshot sums.
- AC2 (checkpoint fine grid, format unchanged): same test reads the final checkpoint (taken at level 2): fine grid, state equal to the last snapshot. The checkpoint code is unchanged; existing checkpoint tests pass.
- AC3 (resume recomputes the level): new `VoxelRun::coarsen_to_fit` (driver.rs; `try_coarsen` now returns bool) is called in `run_with` after `set_coarsening` on resume. Test `coarsened_run_resumes_from_a_mid_run_checkpoint` checkpoints at 1e-8 s at level 2, resumes (level history [(1e-8,1),(1e-8,2)]) and matches the uninterrupted run at 3 snapshots within rtol 1e-6 (atol 1e-12 umol/L). Unit test `coarsen_to_fit_finds_the_level_of_a_fine_state` covers uniform, capped, spiked and off cases.
- Reviewer observation (config summary): `VoxelSimulation::config_summary` now stores `max_coarsen_levels` and `coarsen_tol_factor`. `warn_differences` already tolerates missing keys in an old checkpoint (it logs "absent" and goes on), so old checkpoints still resume, with a warning for these two keys. Test `checkpoint_summary_holds_the_coarsening_settings` checks the stored keys and resuming a checkpoint without them.
- Change beyond the listed items (needed for the 1e-6 criterion): a resumed run started with the Rannacher damping pending, which an uninterrupted run only has right after an injection. With it the resumed snapshots differed from the uninterrupted ones by 87 x the rtol-1e-6 tolerance (atol 0); without it by 1.5e-4 x. Added `VoxelRun::set_rannacher_pending` and, in `run_with`, damping is pending on resume only if a seed was injected at exactly the checkpoint time (`ckpt.seeds` has `time == ckpt.t`). Checkpoint format unchanged. This also tightens resume at level 0; existing resume tests still pass.
- Files touched: src/voxel/driver.rs, src/voxel/sim.rs, tests/voxel_checkpoint.rs. src/voxel/checkpoint.rs and output.rs needed no change.

**Independent review (parent):** opus/high reviewer against base 5a06812: no P0-P2; confirmed the Rannacher-on-resume change correct. Two P3 findings, fixed by the parent in one correction pass: (1) skip `coarsen_to_fit` when resuming at an injection time (damping pending), so the damped step runs on the fine grid as in the uninterrupted run; (2) new tests `resume_at_injection_time_matches_uninterrupted` and `coarsened_resume_at_injection_time_matches_uninterrupted` (factor 1e12; verified to fail without fix 1). Re-review (opus/medium): no P0-P2; one optional P3 (empty-history / snapshot-count guards) applied. fmt, release clippy, full cargo test --release clean. Accepted.

## Blockers / Challenges

The resumed run did not match the uninterrupted one at rtol 1e-6 at first: the cause was the Rannacher damping restarting on resume (measured 87 x the rtol-1e-6 tolerance with it, 1.5e-4 x without it). Fixed as described in Review Notes. Remaining gap: carried rsoda step sizes are not stored in checkpoints (as before), so a resumed run is not bit-identical, only within ~1e-10 relative. The test uses `coarsen_tol_factor = 1e5` to make the short example coarsen after injection; the default factor leaves it at level 0 after the seed.

## Commits

- 584767a feat: Fine-grid outputs and level recomputation on resume (ticket 05)

## Time Spent

12m (ticket-start.js to ticket-complete.js)
