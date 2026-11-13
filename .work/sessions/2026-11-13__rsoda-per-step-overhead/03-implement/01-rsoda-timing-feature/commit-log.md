# Ticket 01: rsoda-timing-feature

**Status:** ✅ Done

## Local Test Result

Ran in `${DEV_DIR}/RUST/rsoda` (log: `.scratch/tests/2026-11-13__rsoda-per-step-overhead/01-verify.log`):
`cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release && cargo test --release --features timing && cargo doc --no-deps --all-features` → exit 0. Without the feature: 117 unit tests, analytic 10, golden 1 (golden files present, not skipped) pass. With it: the same, plus `tests/timing.rs::counts_match_the_solver_counters` ok.

## Review Notes

Base: radkin `70ec63b214d9fc41932e8a03f13904e588f06d84`, rsoda `07a075378603ea7e491432e96bf8b63780305cdc` (main, clean). This ticket's code is entirely in rsoda: branch `perf/per-step-overhead`, commit `08a2abf`. radkin has no change for it.
- Feature `timing = []`, off by default. The timers sit at the two stoda call sites (`wm.prepare`, `wm.solve`), inside `#[cfg(feature = "timing")]`, so every LinearSolver (dense, band, sparse, custom) is covered and the default build has no timer code. The call results are bound before the timer stops and then handled as before (`out?`, `match solved`), so the control flow is unchanged.
- `TimingStats { n_prepare, t_prepare, n_solve, t_solve }` lives in `State` (reset by `start`, i.e. per integration or re-init), and `Lsoda::timing()` reads it. Documented: `t_prepare` includes the Jacobian evaluation (user or FD), so the radkin bench has to subtract its own jac timer.
- The test asserts `n_prepare == nje` on Robertson for the DenseFd and BandedFd kinds. CI gets a `--features sparse,timing` matrix entry (clippy, test, doc). README has a Profiling section.
- rsoda branch pushed (`origin/perf/per-step-overhead`), and radkin's `Cargo.toml` pins rev `08a2abf` (no `[patch]`). radkin `cargo fmt --check`, `clippy -D warnings` and `cargo test --release` pass (54 tests, log `01-radkin.log`).
- Inline ticket, self-reviewed from the diff (localized, additive, cfg-gated).

## Blockers / Challenges

None. Note: `prepare` includes the Jacobian evaluation, so `t_prepare` is not a pure factor time. The bench breakdown (ticket 2) gets the factor share as `t_prepare - jac time`.

## Commits

- f07a9dd feat(deps): pin rsoda 08a2abf with the optional timing feature (ticket 01)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
