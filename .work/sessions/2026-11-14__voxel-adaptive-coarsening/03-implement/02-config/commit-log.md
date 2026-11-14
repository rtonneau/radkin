# Ticket 02: config

**Status:** ✅ Done

## Local Test Result

Command: `cargo test --release voxel::config` (output in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/02-test.txt`).
Result: 7 passed, 0 failed, 54 filtered out (lib tests). Includes `coarsening_defaults`, `parses_coarsening_keys`, `coarsening_divisibility_error_names_axis`, `rejects_invalid_coarsening`, plus the three pre-existing config tests.

Also run:
- `cargo fmt --check`: clean (exit 0).
- `cargo clippy --release --all-targets -- -D warnings`: clean (exit 0, output in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/02-clippy.txt`).

## Review Notes

Base commit: 333d411d4be156af93acb8680ad5cabd312d4e56 (clean tree at start of ticket work, apart from the parent's pre-existing `.session-config.json` change and this untracked commit-log directory).

Acceptance Criteria:
- `VoxelConfig` parses `max_coarsen_levels` (integer >= 0, default 0) and `coarsen_tol_factor` (float > 0, default 0.1) from `[voxel]`: done. Negative, non-integer, string, and out-of-range values (`2^levels` overflows `usize`) are rejected. `coarsen_tol_factor` accepts an integer literal, as the other float keys do, and rejects 0, negatives, and non-finite values.
- Divisibility: done. Checked right after the grid is built, in `VoxelConfig::from_table`, for x, y, z in turn. The error is `RadkinError::Config` with text like `voxel.grid.dims along y (6) must be divisible by 2^voxel.max_coarsen_levels = 4`.
- Unit tests: done, in `src/voxel/config.rs`.

Notes on scope:
- Only `src/voxel/config.rs` was changed.
- The config summary that checkpoints compare on resume lives in `src/voxel/sim.rs` (`RunSettings::summary` and `VoxelRun::config_summary`). It does not cover `[voxel]` grid or diffusion keys, so I did not add the new keys there. The ticket allows only `config.rs`. The parent should decide whether the new keys belong in that summary, probably as a separate ticket. Also, `from_table` does not check `[voxel]` against the summary, so an old checkpoint would not warn about a changed `max_coarsen_levels`.
- `cargo fmt` was run once and reformatted `config.rs` (the new test lines were long). No other files were touched.

## Blockers / Challenges

None

## Commits

- a07a91c feat: Coarsening config keys for the voxel grid (ticket 02)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
