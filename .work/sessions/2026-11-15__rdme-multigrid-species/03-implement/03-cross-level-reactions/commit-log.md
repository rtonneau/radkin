# Ticket 03: cross-level-reactions

**Status:** ✅ Done

## Local Test Result

Before any code: the mesh-correction session had merged (`origin/main` 7dcbccb, PR #16). Rebased `feat/rdme-multigrid-species` onto it (conflicts in `src/ssa/voxel_sim.rs` and `.work/sessions/.current-session` resolved: main's `resolve_correction`, critical-spacing warnings and 5-argument `RdmeSetup::new` kept, the hierarchy built to `hierarchy_depth`; multigrid unit tests use `correct = false`). Full checks after the rebase: exit 0 (`03-rebase-verify.log`).

Verification Step `cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log `.scratch/tests/2026-11-15__rdme-multigrid-species/03-verify.log`): exit 0, 24 test binaries ok. New `tests/ssa_rdme_multigrid.rs`:
- `ssa_multigrid_conversion_to_a_finer_level_matches_the_mean_field` ok: worst deviation 2.23 SEM over the 4 A and 32 B means (400 runs).
- `ssa_multigrid_cross_level_association_has_the_coarse_voxel_rate` ok: mean 4.5345e-9 s vs 1/k_meso 4.5492e-9 s (0.20 SEM, 4000 runs); the uncorrected 1/c (6.2e-9 s) is more than 10 SEM away.

## Review Notes

- Base (after the rebase): `798a423e3f2c4a69f382eb774cdafc3ce3a00c6e`. Inline mode, diff inspected against the criteria.
- The firing rules were written in ticket 02 (`fire_reaction`, `pick_inside`, `uniform_inside`, the coarse-sum updates in `fire_hop`); this ticket removes the cross-level refusal from `run_rdme` and validates them. Cross-level propensities come from the reaction level's `SsaModel`, which carries the ADR 0004 correction for that level's spacing (`RdmeSetup::new` with `correct`).
- Conversion test: A (H_r) on level 1 -> B (H2) on level 0 on an 8x2x2 grid; reference is the mean of the mgRDME generator (Erban & Winkelmann 2024, eq. 21-22, the time-dependent form of eq. 24) integrated with RK4 (dt 1e-10 s), 3D zero-flux Laplacian on both grids, B produced at k/8 into each child of A's voxel.
- Association test: OH_r on level 2 of a 4x4x4 grid (a single 8 nm voxel), e_aq hopping on the 2 nm grid: the waiting time is exactly exponential with the level-2 corrected propensity (hRDME III.C.4, eq. 10-12), independent of the hops; it also checks that the product appears in the electron's voxel and that the hop updates of the coarse rate do not drift (> 10 hops per run).
- `Multigrid::step` made public for the association test (event-by-event loop).
- `MultigridSetup::cross_level_reaction` kept (used by the unit tests, useful to report).

## Blockers / Challenges

The first run of the conversion test failed on a sanity bound of mine (every fine voxel above 0.5 molecules at 1 µs), not on a statistical check; lowered to 0.05 (the far voxels hold 0.09).

## Commits

- e1c9cbc feat: cross-level reactions in the multi-grid RDME, validated (ticket 03)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
