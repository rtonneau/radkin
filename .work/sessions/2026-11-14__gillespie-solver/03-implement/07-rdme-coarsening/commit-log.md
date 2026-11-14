# Ticket 07: rdme-coarsening

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: all `ssa` tests pass, including the 7 unit tests of `ssa::rdme_coarsen` and `ssa_rdme_coarsening_matches_the_fine_grid`.

`cargo test --release --test ssa_rdme ssa_rdme_coarsening -- --nocapture`:

```
coarsening: worst total difference 1.40 combined SEM; hops per trajectory 300906 fine, 59357 coarsened; reactions 7191 and 7209; levels of trajectory 0: [(5.36e-6, 1), (1.47e-5, 2), (1.47e-5, 3), (5e-5, 2), (5e-5, 1), (5e-5, 0), (5.44e-5, 1), (6.58e-5, 2), (6.58e-5, 3)]
test result: ok. 1 passed; 0 failed
```

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 195 passed, 0 failed.

## Review Notes

Base: 4da39187fdd03f01f5e02e1523b4f7633cc4a99a. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- Dispersion test (`Rdme::is_uniform_in_blocks`): aggregated chi-square index of dispersion per species, threshold `z = 30 * coarsen_tol_factor`. Added beyond the ticket text: a per-block test for blocks with at least 16 molecules (threshold `z + sqrt(2 ln B16)`), because the aggregate alone dilutes one cluster in a large uniform background; a unit test shows the cluster is missed without it and caught with it.
- Merge sums the 8 counts; the level's own grid and `SsaModel` (voxel volume) give hop rates and propensities.
- `inject` refines to level 0 first (each molecule to a random child); the integration test uses two seeds and checks the history (level 0 at the second seed time, coarsened again after).
- Conservation: unit test over coarsen x2 and refine x2; the integration test checks that the fine-grid snapshots hold every molecule.
- Snapshots on the fine grid (`fine_field`, in units of 1/8^max_levels molecule so that the ensemble sums stay integers).
- Coarsened vs fine ensemble: 1.4 combined SEM on the species totals, 5 times fewer hops, with the default factor.
- Threshold mapping tuned on the homogenising case and recorded in `docs/ssa.md`: factors 0.03 to 0.3 are indistinguishable from the fine grid (1.3 to 2.6 SEM, 3.8x to 10.6x fewer hops); 1 and 3 are biased (6.3 and 24 SEM). The mapping `z = 30 * factor` with the default 0.1 is kept.

The test runs every 8 x n_voxels events; no test is run at t_start, so a uniform initial solution is not coarsened before the first seed (few events happen before it).

## Blockers / Challenges

None. My first test expected a refinement at the first seed; the grid is still at level 0 there, so the test uses a second seed to exercise it.

## Commits

- cdceb9c feat: RDME coarsening by a Poisson dispersion test, with multinomial refinement (ticket 07)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
