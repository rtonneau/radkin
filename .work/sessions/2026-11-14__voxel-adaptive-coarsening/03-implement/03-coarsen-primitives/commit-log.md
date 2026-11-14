# Ticket 03: coarsen-primitives

**Status:** ✅ Done

## Local Test Result

Verification Step: `cargo test --release voxel::coarsen` -> exit 0, `test result: ok. 12 passed; 0 failed` (lib; other targets 0 tests, filtered). Tests: coarse_grid_halves_dims_doubles_spacing, restrict_conserves_totals, restrict_is_block_mean, prolong_of_restrict_is_exact_for_blockwise_constant, restrict_of_prolong_is_identity, in_place_buffers_are_overwritten, uniform_field_can_coarsen, small_noise_below_tolerance_can_coarsen, single_spike_cannot_coarsen, face_neighbour_blocks_coarsening, single_block_grid_skips_walls, non_finite_cannot_coarsen.

Also run:
- `cargo fmt --check` -> exit 0.
- `cargo clippy --release --all-targets -- -D warnings` -> exit 0, no warnings.
- `cargo test --release` (full suite) -> exit 0, lib 73 passed, every other target ok, no failures.

Output saved in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/evidence/` (`03-test.txt`, `03-clippy.txt`, `03-test-all.txt`).

## Review Notes

Base commit: 8ec0db71e7fae8c939ac488aec16290002bed1e5 (clean tree apart from the session files).

Checked against the Acceptance Criteria:
- AC1: new `src/voxel/coarsen.rs` (registered in `src/voxel/mod.rs`) with `coarse_grid` (dims/2, spacing*2, same origin; `RadkinError::Value` naming the axis on an odd dim), `restrict` (mean of each 2x2x2 block, every species), `prolong` (piecewise constant) and `can_coarsen(grid, state, atol, rtol, factor)` (per species, per block: max - min over the 8 children plus their in-domain face neighbours, up to 24, must be `<= factor * (atol + rtol*|mean of children|)`). Boundary face neighbours are skipped. Both state functions take the fine grid.
- AC2: conservation (3 species, anisotropic spacing, totals = sum * voxel volume, rel. diff <= 1e-14), prolong(restrict(x)) == x exactly for blockwise-constant x, restrict(prolong(y)) == y exactly, uniform field coarsens, single spike (interior, corner and last voxel) does not. Extra: face-neighbour jump between internally uniform blocks blocks coarsening, noise below tolerance passes, NaN/inf fail.
- AC3: `restrict`, `prolong` and `can_coarsen` work on caller buffers and do not allocate. The `restricted`/`prolonged` wrappers allocate the output state once per call.

Design choices for the parent:
- The block sum is pairwise, `((a+b)+(c+d))+((e+f)+(g+h))`, so restricting eight equal values is exact. That is what makes restrict(prolong(y)) == y bit-exact.
- The test is `spread <= tol` rather than strict `<`, so a field of zeros coarsens even with atol = 0. Non-finite values fail through the block mean, since every fine voxel is a child of some block.
- Size and odd-dim mismatches in `restrict`/`prolong`/`can_coarsen` panic (`assert!`, once per call). They are programmer errors: config (ticket 02) already rejects dims not divisible by 2^L.
- `src/voxel/grid.rs` is unchanged. Nothing needed to be added there. No re-exports were added to `mod.rs`. No driver wiring (ticket 04).

Review rounds: 0 (no external review in this pass).

## Blockers / Challenges

None. One test threshold was borderline on the first draft (the spike in `small_noise_below_tolerance_can_coarsen` gave a spread equal to the tolerance). I raised it to 2e-2 before the first run.

## Commits

- 929f08c feat: Restrict, prolong and coarsening test for voxel grids (ticket 03)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
