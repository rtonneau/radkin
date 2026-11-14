# Ticket 01: unit-cell-exact-ssa

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: 9 unit tests (`ssa::rng`, `ssa::model`) and 5 integration tests (`tests/ssa_exact.rs`) pass, 0 failed.

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 148 passed, 0 failed.

## Review Notes

Base: 401fc51028a9f288f3ae58aa104caa6c054ac76c. Pre-existing changes: none outside `.work/` (session files owned by gps).

Checked against the Acceptance Criteria:
- Propensity model (`src/ssa/model.rs`): mass action with falling factorials, fixed species as concentration factors, Michaelis-Menten rate law, `need` list for every species consumed on balance (negative products included).
- A + A: the stored rate constant already carries radiopyo's convention (a scalar `k_value` of a single-reactant reaction is divided by its stoichiometry in `apply_nk_convention`), so the model takes `mat.k.k_values` as is: `a = k/(N_A V) n (n-1)`, two molecules removed. Unit test for `2 OH_r -> H2O2` and `H_r + H_r -> H2`; the integration test shows a factor-2 error would be more than 50 SEM away.
- RNG: ChaCha8, stream = trajectory index; uniform, exponential and integer samplers use only IEEE basic operations (`det_ln` instead of libm).
- Exact trajectory: direct method, output rows at `t_eval`, identical for the same seed and stream.

Changed after the first run: integration test names got an `ssa_` prefix, because the Verification Step filter `cargo test --release ssa` matches test names and skipped them otherwise. The dose-rate segment and source handling are already in the stepper core; the beam schedule and its tests are ticket 2.

Only `rand_chacha` is added (it re-exports `rand_core`); the `rand` crate is not needed because the samplers are in `src/ssa/rng.rs`.

## Blockers / Challenges

None. Two test expectations of mine were wrong at first (H2O is a built-in fixed species; the stored k of `2 A` is `k_value / 2`) and were corrected against the code.

## Commits

- dca3a07 feat: unit-cell exact SSA with a propensity model and deterministic RNG (ticket 01)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
