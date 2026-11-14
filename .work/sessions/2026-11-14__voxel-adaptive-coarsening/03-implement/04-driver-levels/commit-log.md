# Ticket 04: driver-levels

**Status:** ✅ Done

## Local Test Result

Verification Step: `cargo test --release voxel::` -> exit 0, lib `test result: ok. 57 passed; 0 failed` (the other targets have no test that matches the filter). New tests: `homogenising_run_coarsens_and_keeps_the_totals`, `injection_after_coarsening_refines_and_conserves_the_deposit`, `set_coarsening_checks_the_grid_and_refines_above_the_maximum`. All 7 existing driver tests pass unmodified.

Also run:
- `cargo fmt --check` -> exit 0.
- `cargo clippy --release --all-targets -- -D warnings` -> exit 0, no warnings.
- `cargo test --release` (full suite) -> exit 0. Lib 76 passed (73 before, plus the 3 new tests). Every integration target is ok (voxel_checkpoint 3, voxel_geant4_seed 3, voxel_seed 1, golden_runs 5, ...), no failures.
- Mutation check: with `fine_per_voxel` forced to 1.0 (totals not weighted by the level's volume), both new conservation tests fail. The file was restored afterwards.
- CLI smoke test (scratch copies of `data/voxel_example.toml` run to `t_end = 1e-5`, coarsening off and on with `max_coarsen_levels = 2`). Both exit 0, 58 steps each. The "on" run logs `coarsening: 5 level changes, final level 2`: 0->1->2 on the uniform start, refined to 0 by the seed at 1e-12 s, then 0->1->2 once homogenised. Final `yields.csv` molecule totals agree within 5.5e-10 relative (max over species). Snapshots stay on the fine grid (768 values per field). Voxel half-steps drop from 222720 to 193320. This is a single run, not a timing measurement; ticket 06 does that.

Output saved in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/evidence/`: `04-fmt.txt`, `04-clippy.txt`, `04-test-voxel.txt`, `04-test-all.txt`, `04-cli-off.txt`, `04-cli-on.txt`, `04-compare.txt`.

## Review Notes

Base commit: 19d13b10b76eeb8cf8a9e2580a43a051633878be (clean tree apart from the session files).

Checked against the Acceptance Criteria:
- AC1: `VoxelRun` has a private `level`, plus the grid of every level up to the maximum (`grids[0]` is the fine grid). `state`, `diffusion` and `carry` live on `grids[level]`. Also new: `level()`, `fine_grid()`, `fine_state()` and the pub `level_history: Vec<(f64, u32)>` (`(t, level)` after each change).
- AC2: in `advance_to`, after each accepted step (after the dt update), `try_coarsen` runs. If `level < max_coarsen_levels` and `can_coarsen(grids[level], state, control.atol, control.rtol, coarsen_tol_factor)` passes, the state is restricted, `Diffusion` is rebuilt for the coarse grid with the same coefficients, and `carry.reset_voxels(n)` sizes the per-voxel steps for the new count and zeroes them, keeping the background step. It drops one level per accepted step. dt is not touched.
- AC3: `inject` calls `refine_to(0)` (piecewise-constant prolongation, level by level) before the existing assert and add. The Rannacher flag, the dt reset and the carry reset are unchanged. The test checks `rannacher_pending`, `dt() == dt_initial`, that the deposit sits exactly on its fine voxels on top of the prolonged state, and that B + C conserves the deposit (1e-9) after reaction and diffusion.
- AC4: with the default `max_coarsen_levels = 0`, `try_coarsen` returns at once and `domain_sums` multiplies by `8^0 = 1.0` (exact), so default runs are bit-identical. No existing test was changed.
- AC5: the homogenising test (A + B -> C, 8x4x4 grid, t_end 3e-6 s) reaches level 2 at about 1.7e-6 s and keeps taking steps there. Domain totals match the uncoarsened run within rtol 1e-6 (measured about 2.5e-10). A + C and B + C hold at every totals record through the level changes (1e-9). The fine field stays within the control tolerance of the fine run (measured scaled error 4e-5). The coarsened run does fewer voxel half-steps. The injection test is described under AC3.

Notes items:
- Domain totals keep their unit (µmol/L summed over fine voxels): `domain_sums` and the background totals weight each voxel by `8^level`. That is the current voxel volume over the fine voxel volume (exact in floating point), so `yields.csv` and the G-values need no change.
- The background vector is grid-independent and carries over. A uniform field restricts to exactly the same values (pairwise mean), so background voxels stay background on every level. The test checks this: `background_voxel_steps == voxel_steps` at level 2.
- `voxel_steps` and `background_voxel_steps` count the voxels of the level they ran on, i.e. the work done.

Interim for ticket 05 (parent note: no caller should silently get a coarse state):
- `sim.rs` now calls `run.set_coarsening(voxel.max_coarsen_levels, voxel.coarsen_tol_factor)`, so the `[voxel]` keys take effect. Snapshots and checkpoints read `r.fine_state().into_owned()` instead of `r.state.clone()`. Outputs and checkpoints are therefore on the fine grid already, and the checkpoint format is unchanged. `VoxelOutput` gains `level_history` (the run stats for ticket 06). The two `VoxelOutput` literals in the `output.rs` tests needed that field to compile; that file is not in this ticket's list, and it is a 2-line change.
- Left for ticket 05: resume does not yet recompute the level from the checkpoint state. A resumed run starts at level 0 and coarsens again one level per accepted step. Ticket 05 also owns the snapshot, output and checkpoint tests at level > 0. `run.state` stays pub and holds the current level's state; the only non-test readers (snapshots, checkpoints) now use `fine_state()`.
- `set_coarsening` validates the factor (finite, > 0) and that the fine dims are divisible by `2^max_levels` (naming the axis), before it changes anything. Lowering the maximum below the current level refines to it.

Allocation: restrict/prolong, the coarse grids and the rebuilt `Diffusion` allocate only at level changes, `fine_state()` only at output time. `can_coarsen` does not allocate. The kernel is untouched.

Review rounds: 0 (no external review in this pass).

**Independent review (parent):** fresh opus/xhigh reviewer against base 19d13b1, stable diff: no P0-P3 findings. Observations: (1) a step landing on a seed time may coarsen then immediately refine in inject (harmless, noted for ticket 06 level counts); (2) checkpoint config_summary omits the coarsening keys (carried to ticket 05). Gaps closed by the parent: debug clippy and `--features timing` clippy both clean. Accepted: all ACs met, credible verification, no unresolved actionable findings.

## Blockers / Challenges

None. Two small points on the way:
- A borrow conflict in `record_totals`: `fine_per_voxel()` was called while `background_totals` was borrowed mutably. I moved the computation before the borrow.
- In the first calibration only 1-2 steps ran after the grid reached level 2, because dt grew too fast. The test helper sets `dt_max = 1e-7` so that the coarse levels are actually exercised (about 13 steps on level 2).

## Commits

- c369f5e feat: Global coarsening levels in the voxel driver (ticket 04)

## Time Spent

19m (ticket-start.js to ticket-complete.js)
