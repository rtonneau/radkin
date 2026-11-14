# Ticket 02: species-level-state-and-hops

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log: `.scratch/tests/2026-11-15__rdme-multigrid-species/02-verify.log`): exit 0, every binary ok (lib: 129 passed). New: `ssa::rdme_multigrid::tests::{levels_of_reactions_follow_the_coarsest_reactant, a_coarse_species_diffuses_on_its_level, a_coarse_reaction_has_the_macroscopic_rate, all_fine_levels_match_the_single_level_rdme, rates_and_level_sums_stay_consistent_across_levels, poisson_fill_and_fine_field_use_the_species_level, injection_goes_to_the_containing_voxel}`; `ssa_rdme_species_levels_are_checked` extended. Existing fixed-seed RDME tests (`tests/ssa_rdme.rs`, coarsening and tau) pass unchanged.

## Review Notes

- Base: `7ada7b0fa34190453da2bd957acd08b465efe724`. Inline mode, diff inspected against the criteria (no independent reviewer: the harness avoids spawning agents unless asked).
- New `src/ssa/rdme_multigrid.rs`: `MultigridSetup` (species and reaction levels, per-level readers, hoppers, parent maps, site offsets) and the `Multigrid` trajectory. Counts of each species live on its level, with sums at every coarser level up to the highest level whose reactions read it; one next-subvolume queue over the sites of every level. Hops update the coarse sums only when they cross a parent boundary (incremental reader update, as in `Rdme`); reaction events recompute the touched sites in full.
- The cross-level firing rules (pick a contained voxel in proportion to its counts; product at the finest consumed molecule's spot, uniform below it) are already in `fire_reaction`, because the event structure needs them; runs with a cross-level reaction are still refused in `run_rdme` until ticket 03, as planned. `rates_and_level_sums_stay_consistent_across_levels` exercises them for bookkeeping only (sums, rates, conservation), not their statistics.
- Product placement refinement against the grill wording: a product goes to the spot of the finest consumed molecule (uniform below it if the product is finer), which reduces to "uniform among the children" (mgRDME eq. 20) when the reactant sits on the reaction level, and keeps e.g. `e_aq + O2(coarse) -> O2_minus` at the electron's voxel. To be stated in the docs (ticket 06).
- Coarse diffusion test: exact mean of the 1D zero-flux hop process (cosine eigenmodes) on the level-1 4x1x1 grid, 200 runs, 4 SEM.
- `run_rdme`: dispatch via `enum Model { Single, Multi }` and a private `VoxelTrajectory` trait; species levels all 0 -> the old `Rdme` path, same RNG use. Hierarchy built to `hierarchy_depth`; the field unit uses `8^depth`. Temporary config errors: species levels with `max_coarsen_levels > 0` (ticket 04), with `tau_leap` (ticket 05), and cross-level reactions (ticket 03).
- `RdmeOutput::species_levels` and `snapshots.json` key `species_levels` (species order).
- `rdme.rs`: `strides` and `face_weight` made `pub(crate)`.

## Blockers / Challenges

None. Next ticket (03) requires the rebase onto `main` after `2026-11-15__rdme-mesh-corrected-rates` merges.

## Commits

- f93cfc4 feat: multi-grid RDME state, hops and same-level reactions (ticket 02)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
