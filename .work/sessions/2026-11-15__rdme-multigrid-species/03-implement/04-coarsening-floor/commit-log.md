# Ticket 04: coarsening-floor

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log `.scratch/tests/2026-11-15__rdme-multigrid-species/04-verify.log`): exit 0, 24 test binaries ok. New: `ssa::rdme_multigrid::tests::{global_levels_are_a_floor_and_conserve_every_species, the_dispersion_test_reads_only_the_species_that_merge, a_mixing_grid_coarsens_and_a_seed_refines_it}`, `ssa_rdme_species_levels_follow_the_global_coarsening` (H2 on level 1, global coarsening to 3, 24 trajectories: worst 1.58 SEM against the fine grid, 1.26e6 hops vs 7.05e6 fine). `ssa_rdme_coarsening_matches_the_fine_grid` (single-level path, refactored dispersion test) unchanged and ok.

## Review Notes

- Base: `7aba38b7dc3815fe192791d37c7862c12342c3ca`. Inline mode, diff inspected against the criteria.
- `MultigridSetup` now holds one `Layout` per global level `g` in `0..=max_global` (`voxel.max_coarsen_levels`): species level `max(g, own)`, reaction levels, readers and hoppers recomputed per layout. `Multigrid` holds `g` and switches layout in `set_global`: species whose level rises are summed into their parents, those whose level falls are split molecule by molecule into a random child (one level at a time), then everything is rebuilt.
- Dispersion test: `rdme_coarsen::blocks_are_uniform` extracted as a pure function (same arithmetic and order); `Rdme::is_uniform_in_blocks` calls it with every species, `Multigrid::is_uniform_in_blocks` with the species living on level `g` only (those that would merge). Cadence `8 x voxels of level g` events, as in `Rdme`.
- `inject` refines to `g = 0` first (each species back on its own level). `level_history` reported in `RdmeOutput` for multi-grid runs too. The `max_coarsen_levels > 0` refusal in `run_rdme` is removed.
- Finding while writing the end-to-end test: putting OH- (a track product) on level 1 deviates strongly from the fine grid (+40 % H3O+ at 1 us, 30+ SEM), with or without global coarsening, so it is not the floor. The seeded cluster straddles the 200 nm voxel boundary, OH- is spread over about 8 times the volume of its fine H3O+ partner, and their near diffusion-limited recombination slows down. That is the resolution trade-off of ADR 0005, not a defect (the conversion and association tests validate the rules); the test uses H2 and its comment records the OH- case. To be documented in `docs/ssa.md` (ticket 06).

## Blockers / Challenges

One wrong expectation of mine in a unit test (level-1 voxel count of an 8x4x4 grid is 16, not 32), fixed.

## Commits

- 984929b feat: global coarsening as a floor for species levels (ticket 04)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
