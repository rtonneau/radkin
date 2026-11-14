# Ticket 01: species-levels-config

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log: `.scratch/tests/2026-11-15__rdme-multigrid-species/01-verify.log`): exit 0, every test binary `ok` (lib: 122 passed). New tests ok: `ssa::config::tests::{parses_species_levels, species_levels_need_a_stochastic_voxel_run, depth_covers_coarsening_and_species, suggestion_follows_the_square_root_of_the_diffusion_ratio}`, `ssa_rdme_species_levels_are_checked`.

## Review Notes

- Base: `c4455726d12c57485757fd964674c7b44288e015` (plan commit). Pre-existing changes: only gps session files (`.session-config.json`, `03-implement/`), owned by the scripts.
- Inline mode (no independent reviewer: the harness avoids spawning agents unless asked); diff inspected against the criteria.
- Parsing: `[ssa.species_levels]` -> `SsaConfig::species_levels` (label, level) in file order; non-table, negative, non-integer, string and > 16 values are errors naming `ssa.species_levels[.<name>]`. Unknown species are caught when resolved against the dynamic species (`resolve_species_levels`), in `run_rdme`.
- Refusals: `require_volume` (unit-cell SSA), `run()` in the CLI (any unit-cell run), and `voxel()` for the deterministic path (`refuse_species_levels`).
- Depth: `hierarchy_depth`; the divisibility check is now `voxel::config::check_divisible(dims, levels, key)`, shared with `voxel.max_coarsen_levels` (same message format), and names `ssa.species_levels.<deepest species>` when species levels set the depth.
- Suggestion: pure `suggested_level(d, d_slowest, depth)`; `run_rdme` logs it (info) for species with a non-zero initial concentration that are not seeded, capped by the deepest level the grid dimensions divide into. Seeded species (e.g. H3O+ from the track) are skipped: they carry the track structure.
- Until ticket 02, a level > 0 is a config error ("not implemented yet") so that a run never silently ignores it; level 0 everywhere runs the plain RDME. Ticket 02 removes that refusal.
- `SsaConfig` literals in `tests/ssa_tau.rs` and `examples/ssa_bench.rs` got `..Default::default()`.

## Blockers / Challenges

None.

## Commits

- 0286a80 feat: [ssa] species_levels table, depth and level suggestion (ticket 01)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
