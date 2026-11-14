# Ticket 01: rsoda-restart

**Status:** In Progress

## Local Test Result

Ran in `$DEV_DIR/RUST/rsoda` (`D:\DEV\RUST\rsoda`, branch `feat/lsoda-restart` from `f5d1312`, changes uncommitted):

```bash
cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features restart
```

Result: exit 0. fmt and clippy silent. Filtered tests: 3 lib tests (`lsoda_restart_keeps_state_storage_and_linear_solver`,
`lsoda_restart_reallocates_only_for_another_storage`, plus the existing `stoda_repeated_error_test_failures_restart_at_order_one`)
and 5 in `tests/restart.rs` (`restart_matches_new_dense_and_banded`, `restart_matches_new_sparse`,
`restart_keeps_a_custom_linear_solver`, `restart_after_an_error_and_a_jacobian_change`,
`restart_rejects_bad_input_without_change`) all pass. Output: `.scratch/tests/2026-11-14__rsoda-lsoda-restart/evidence/01-verification.txt`.

Extra checks (`evidence/02-full-suite.txt`, exit 0): `cargo test --release --all-features` (all suites, 122 lib tests incl. the new
`stoda_state_reset_matches_a_new_state_and_keeps_storage`), `cargo clippy --release --all-targets -- -D warnings` and
`cargo test --release` without features (sparse tests compiled out), `cargo doc --no-deps --all-features` without warnings.

## Review Notes

Checked against the Acceptance Criteria:

- AC1: `Lsoda::restart(&mut self, t0, y0, h0) -> Result<(), LsodaError>` added (`src/lsoda.rs`). It validates `y0.len()` and
  `h0` (NaN), sets `opts.h0 = h0`, `t`, `y` (copied into the existing vec), clears `init`/`changed`/`halted`/`nhnil`, and resets the
  state in place. `State::reset` (`src/state.rs`) rebuilds the state through the same private `State::build` that
  `State::with_solver` uses, with the work arrays (`y`, `yh`, `ewt`, `savf`, `acor`) moved out, cleared and zero-resized in their
  own storage and the linear solver moved back in (a zero-sized `Detached` placeholder, no allocation, holds its slot). Sharing
  the constructor means a field `new()` sets cannot be forgotten by the reset. Block C (`start`) now uses the same in-place reset,
  so the first call after `new` no longer reallocates the state and solver `create` already built either. The solver is replaced
  only when its storage no longer fits `opts.jacobian` (after `set_jacobian`), through the existing take/release pool (custom
  solver kept). Unit tests check by pointer that `yh` and the boxed linear solver survive `new` -> run -> `restart` -> run for the
  dense and banded kinds, and that only a storage change swaps the solver.
- AC2: `tests/restart.rs` runs one solver through four Robertson states `(t0, y0, h0)` (computed and given `h0`, t0 from 0 to 4e3),
  each integrated over ten decades so BDF and Jacobians are used, restarting between states, and compares the bits of every output
  `t`/`y` plus `nst`, `nfe`, `nje`, `nqu`, `nqcur`, `hu`, `tcur`, `tsw`, `mused` with a fresh `Lsoda::new` (same options, `h0`).
  Identical for DenseUser, DenseFd, BandedUser, BandedFd, SparseUser, SparseFd, a custom `SparseLu` via `with_linear_solver`, and
  after a `MaxStepsExceeded` error plus a `set_jacobian` to another storage. The state right after `restart` also reads as new
  (t, y, tcur, tsw, counters, method, order).
- AC3: fmt, clippy `-D warnings` (with and without features) and tests clean. Commit and push are left to the parent.

Side effect worth noting for review: the old `start` swapped in a placeholder solver and only stored the new state on success, so a
start that failed after taking a custom solver (e.g. RHS error at `t0`) dropped the user's solver and the next start would panic
in `take_solver`. The in-place reset removes that path. On a failed start the state is now partly reset (e.g. `nfe = 1` after the
initial RHS call) instead of left as before; `init` stays false, so the next call starts again from `t`/`y` as before, and all
existing tests pass. Restart assumes the system's sparsity pattern is unchanged (documented on `restart`), as `set_jacobian`
already does.

Review rounds: 0 (self-review only).


### Parent review record

- Base: radkin a91c033300abb3ca0b4281536ef161f0043cd5f9 (no radkin code changes); rsoda base f5d1312, work on rsoda branch feat/lsoda-restart.
- Independent reviewer (fresh context, opus/high): REVIEWED, 2 findings, no P0-P2.
  - P3 src/lsoda.rs:2050 sparse kinds lack a storage-pointer test. Disposition: accepted without change; the sparse solver goes through the same same_storage/move-back path the reviewer verified, the bitwise sparse tests pass, and the rober test helper has no sparsity pattern. Ticket 03's bench will show any per-voxel reallocation.
  - P3 src/lsoda.rs:292 restart doc overstated no reallocation. Fixed inline by the parent.
- Gap: radkin golden tests vs the new rsoda not run yet; covered by ticket 03 (rev bump runs the full radkin suite).
- Parent re-ran rsoda fmt, clippy -D warnings --all-features, full tests --all-features: clean. Committed and pushed rsoda 5c8c501 on origin/feat/lsoda-restart.
- Accepted: all criteria met.

## Blockers / Challenges

None. One test of my own needed changing: an accuracy check of a banded FD run integrated straight to 4e5 against a
tight-tolerance reference failed the 1e-4 relative bound (a property of that direct run, not of restart); it now compares bitwise
with a fresh banded FD solver instead.

## Commits

Filled by ticket-complete.js.

## Time Spent

Filled by ticket-complete.js.
