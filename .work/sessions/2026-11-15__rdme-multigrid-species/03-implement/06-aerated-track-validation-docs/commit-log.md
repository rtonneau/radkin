# Ticket 06: aerated-track-validation-docs

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log `.scratch/tests/2026-11-15__rdme-multigrid-species/06-verify.log`): exit 0, 25 test binaries ok. New `tests/ssa_rdme_aerated.rs::ssa_aerated_track_with_coarse_o2_matches_the_fine_grid` (16 trajectories, fixed seed, 0.3 s): O2 on level 1 worst 2.95 SEM (4.1x fewer hops), level 2 worst 3.74 SEM (14.2x fewer hops), all within 4 combined SEM + 1 molecule.

Measurement for the docs (`cargo test --release --test ssa_rdme_aerated -- --ignored --nocapture`, 96 trajectories): sparse track level 0/1/2: 200 270 / 48 584 / 14 121 hops, 0.93 / 0.30 / 0.14 s, worst 3.42 and 2.22 SEM; dense track: 205 682 / 54 642 / 20 412 hops, worst 22.0 and 28.5 SEM.

## Review Notes

- Base: `553e400d104f266f6d0474016175c082b12922fa`. Inline mode, diff inspected against the criteria.
- First version of the test used a dense track (1500 each of e_aq, OH, H3O+, H in a 140 nm cube) and failed at level 1 (e_aq 443 vs 473 at 39 ns). Cause: the track outnumbers the O2 of its voxels (about 1200 O2 in the 8 track voxels vs 6000 radicals); the fine grid depletes O2 locally and a coarse O2 voxel spreads the depletion. This is the accuracy limit the grill named ("a too-coarse O2 level smears O2 depletion"), so the CI test uses a sparse track where O2 is a background, and `docs/ssa.md` reports both regimes. No code change.
- `docs/ssa.md`: new "Species levels" section (key, levels and depth, reaction and product rules incl. the finest-consumed-molecule placement, floor with global coarsening, tau-leaping, suggestion, validation table, aerated measurement table, and the dense-track and OH- counterexamples). README: `[ssa.species_levels]` in the `[ssa]` example and a bullet in the voxel SSA list. `data/voxel_example.toml`: commented `[ssa.species_levels]` example. AGENTS.md layout line for `src/ssa/`.
- Not touched: `examples/voxel_bench.rs` (listed in the ticket): it benchmarks the deterministic voxel path only; the multi-grid measurement is the ignored test `ssa_aerated_track_measurement`, named in the docs.
- Changelog fragment: written by `/gps finish` (`changelog-prepare.js`), as for the earlier sessions.

## Blockers / Challenges

The dense-track failure above (accuracy limit, documented).

## Commits

- 8c085b5 test: aerated-track validation of species levels; docs (ticket 06)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
