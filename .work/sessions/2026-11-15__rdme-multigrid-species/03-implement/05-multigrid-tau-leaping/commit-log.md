# Ticket 05: multigrid-tau-leaping

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` (log `.scratch/tests/2026-11-15__rdme-multigrid-species/05-verify.log`): exit 0, 24 test binaries ok. New `ssa_rdme_multigrid_tau_leaping_matches_the_exact_multigrid` (large-count case, OH_r background on level 1, 16 trajectories, `tau_epsilon = 0.03`, fixed seed): totals within 4 combined SEM + epsilon x peak at every record time, worst difference 0.138 % of the species peak; standard deviations within 80 % + 3 molecules of each other; events per trajectory 3 953 883 exact vs 3 954 854 tau; wall time 3.05 s exact vs 1.43 s tau (2.1x). The single-level tau test is unchanged and ok.

## Review Notes

- Base: `c8a3c073e7cee98abdded0eeb92963f1dc6e2915`. Inline mode, diff inspected against the criteria.
- New `src/ssa/rdme_multigrid_tau.rs` (`MultigridTau`, `impl Multigrid`), a mirror of `rdme_tau.rs`: populations are (species, voxel of its current level); channels are the reactions of every site of `by_level` and the per-level hops; criticality, the CGP bound, Poisson firings, rejection-and-halving, exact batches (threshold `10 + 4 x sites`) and the coarsening count work as in the single-level leap. After a leap the own counts are copied back, the level sums rebuilt and the queue marked stale (rebuilt before the next exact step; `rebuild` clears the flag).
- Same-level reactions change their voxel's counts directly. A cross-level reaction is applied firing by firing with the exact path's rules (`pick_inside` on the start-of-leap counts, products at the finest consumed molecule's spot); its moments in the bound are spread over the contained voxels in the same proportions (`add_reaction_moments`).
- The `tau_leap` refusal in `run_rdme` and its test are removed; `run_rdme` calls `Multigrid::set_tau_leaping`. Fields and helpers of `Multigrid`/`MultigridSetup`/`Layout` that the module needs are now `pub(crate)`.

## Blockers / Challenges

None.

## Commits

- ccb0814 feat: tau-leaping with species levels (ticket 05)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
