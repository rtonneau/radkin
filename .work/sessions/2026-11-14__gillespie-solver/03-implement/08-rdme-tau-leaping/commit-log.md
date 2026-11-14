# Ticket 08: rdme-tau-leaping

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: all `ssa` tests pass, including the 3 unit tests of `ssa::rdme_tau` and `ssa_rdme_tau_leaping_matches_the_exact_rdme`.

`cargo test --release --test ssa_rdme -- --nocapture`:

```
voxel tau-leaping: worst total difference 0.118 % of the species peak; events per trajectory 4103369 exact, 4103020 tau; wall time 2.33 s exact, 1.40 s tau (1.7x)
test result: ok. 6 passed; 0 failed
```

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 199 passed, 0 failed. `cargo test --release --no-default-features ssa`: 66 passed, 0 failed.

## Review Notes

Base: 5434515712efd2f7e964f87c82f0340d26b1a1ee. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- `src/ssa/rdme_tau.rs`: reactions and hops (per species and face) leaped with Poisson counts; a channel is critical when the population it consumes is below 10 and then fires one at a time (`pick_critical`).
- Step: the CGP bound over every (voxel, species) population, hops in and out in the mean and variance (`plan_leap`). The fallback differs from the ticket text on purpose: a leap reads every voxel, so the threshold is `10 + 4 x voxels` expected non-critical events instead of 10/a0, and exact events are taken in batches that double (8 to 512 per voxel) while leaping is refused.
- Non-negative counts (trial state, rejection and halving), molecules conserved by hops: unit test over 10 streams with three conservation laws.
- A leap never crosses the requested time (unit test checks `t` exactly); coarsening counts the leaped events (unit test with both on).
- Tau vs exact RDME ensembles: totals within 4 combined SEM + epsilon x species peak (measured 0.12 % of the peak), mean field of every voxel likewise. As in ticket 5 the criterion carries the epsilon term, because tau-leaping has a systematic error proportional to epsilon.
- Speed-up recorded in `docs/ssa.md`: 1.6x on the seeded test case, 5.1x and 38x on uniform grids with 3000 and 81 000 molecules per voxel.

Changed after the first measurement: on the sparse `data/voxel_example.toml` the first version was 6.4 times slower with `tau_leap = true` than without (one leap decision per output time). The refusal state is now kept across calls (`exact_left`): 1.2 s against 1.0 s, with byte-identical `yields.csv`.

`AGENTS.md` gained a line for `src/ssa/` (not in the ticket's file list).

## Blockers / Challenges

None.

## Commits

- f8d66f7 feat: tau-leaping of reactions and hops in the voxel RDME (ticket 08)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
