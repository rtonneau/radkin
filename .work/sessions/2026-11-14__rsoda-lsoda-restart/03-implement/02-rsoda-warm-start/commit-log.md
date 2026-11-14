# Ticket 02: rsoda-warm-start

**Status:** In Progress

## Local Test Result

Run in `${DEV_DIR}/RUST/rsoda` on `feat/lsoda-restart` (base 5c8c501), changes uncommitted:

```bash
cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features
```

Exit 0. fmt and clippy are clean. Tests: lib 123 passed (one new stoda test), analytic 10, golden 1,
restart 9 (5 from ticket 01, unchanged and passing, plus 4 new), sparse 5, timing 1, doctest 1.
0 failed. Full output: `.scratch/tests/2026-11-14__rsoda-lsoda-restart/evidence/02-verification.txt`.
Also checked: `cargo clippy --release --all-targets -- -D warnings` and `cargo test --release` with
default features are clean, and `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features`
builds.

## Review Notes

**Decision: only the method is carried. The order is not.** `Lsoda::restart_warm(t0, y0, h0)`
starts the new problem on the solver's current method (`mcur`) at order 1.

Why order > 1 is not sound at a restart:
- DSTODA's first step (`Start`) always sets `NQ = 1`. A restart has only `y0` and `h0*f(t0, y0)`,
  so Nordsieck columns 3 and up are zero. At order q > 1 the step is in fact an order-1 theta
  method (`y1 = y0 + h(1 - l0) y0' + l0 h f(y1)`), but the order-q error estimate and step
  exponent `1/(q+1)` judge it. DSTODA's own L640 comment covers this case: after 3 error-test
  failures it assumes the history "has errors of the wrong order" and drops back to order 1.
- Prototype (a temporary order hook in stoda Start plus counters for error-test failures and L640
  rescues, since reverted). Robertson and HIRES, rtol 1e-6, previous solve ending on BDF at q3 or
  q4. There were 28 restart cases: with and without a diffusion-like perturbation of y, and
  `h0 = 0.01*hu` or auto. Compared against a 1e-12 reference
  (`evidence/02-proto-summary.txt`, raw `02-proto-order-vs-method.txt`, code `02-proto.diff` and
  `02-proto_warm_test.rs`):
  - BDF at the carried order: sum nfe 10536, 287 error-test failures, 23 order-1 rescues
    (in 11/28 cases), max err/tol 5.22.
  - BDF at order 1 (method only): sum nfe 9966, 254 failures, 3 rescues, max err/tol 5.00.
  - The carried order needs more nfe than method-only in 20/28 cases. Accuracy is about the same
    because the error control catches the bad history, but it costs more. This is the "error-test
    failures" condition in the AC, so only the method is kept. The `restart_warm` doc comment
    says this and gives the reason.

Acceptance Criteria:
- Opt-in warm start, default unchanged: `restart_warm` is a new method, and `restart` now calls
  `restart_on(Method::Adams, ..)`, which runs the same code as before. The only change in stoda
  (load the BDF coefficients last when `METH = 2` at Start) is skipped on Adams. All 5 bitwise
  `restart` vs `new` tests from ticket 01 and the golden tests pass unchanged.
- Order fallback documented: yes, see above. A test checks `nqcur() == 1` after `restart_warm`
  (the previous order was > 1) and `nqu() == 1` after the first step.
- Stiff test within tolerance of a fresh solve:
  `warm_start_stays_within_tolerance_{dense_and_banded,sparse}` cover all 6 Jacobian kinds. It
  runs Robertson over ten decades from the 4 ticket-01 states. At every output, warm and fresh are
  each within 5x tol of a tight reference (measured max: warm 2.75, fresh 2.09), and warm vs fresh
  is within 10 (measured 4.5; `evidence/02-warm-diffs-measure.txt`). The test also checks that the
  first step is BDF, order 1, with a Jacobian formed and no method switch reported.
- Other new tests: `warm_start_on_adams_and_cold_restart_after_a_warm_one` (on Adams the warm
  start is bitwise the cold one; a `restart` after a warm start is bitwise a new solver) and
  `warm_start_rejects_bad_input_without_change`. The stoda unit test
  `stoda_first_step_on_bdf_leaves_the_bdf_coefficients_loaded` checks the BDF ELCO/TESCO, CM1/CM2,
  LMAX = 6 and NJE = 1.
- Branch pushed: not done (the parent commits and pushes).

Implementation: a `start_method` field on `Lsoda`, set by `restart_on`. `reset_state` applies it:
`METH = BDF`, `MITER = JTYP`, `MAXORD = MXORDS`, the values of a switch to BDF (statement 478 and
block F). This also covers repeated restarts of the problem before the first completed step.
`mcur()` reports BDF right after `restart_warm`.

Suggested commit-message note: "carries only the method; order stays 1 because a restart has no
Nordsieck history (higher columns zero); starting at the previous order gave more error-test
failures, DSTODA order-1 rescues and nfe with no accuracy gain (Robertson, HIRES)".

### Parent review record

- Base: radkin a91c033 (session record commits since); rsoda base 5c8c501, branch feat/lsoda-restart.
- Independent reviewer (fresh context, opus/high): REVIEWED, 1 finding P3, no P0-P2. Confirmed the default restart/new paths are bitwise unchanged (the stoda BDF-start branch is unreachable without restart_warm), BDF-at-order-1 start matches a normal switch, and the order-not-carried claim is supported.
  - P3 prototype summary missed 1 of 28 cases. Fixed by the parent: numbers above corrected from the raw output (28 cases, 10536/287 vs 9966/254, 20/28).
- Gaps accepted: no end-to-end warm BDF->Adams switch-back test (the stoda unit test covers CM1/CM2); HIRES only in the prototype.
- Parent committed and pushed rsoda on origin/feat/lsoda-restart.
- Accepted: all criteria met; the method is carried, the order is not.

## Blockers / Challenges

- `src/stoda.rs` is not in Files to Touch, but it had to change. DSTODA's Start loads the Adams
  coefficients last, so a start on BDF would have used Adams ELCO/TESCO for its order changes.
  This is a 3-line guarded branch plus one unit test. Adams behaviour is unchanged.
- First draft of the stiff test: asserting `nje == 1` after the first step was wrong. With
  `h0 = 0` (auto) the large initial step fails to converge and forms the matrix again. Changed to
  `nje >= 1`.
- For radkin (not part of this ticket): in the prototype, a cold Adams restart with
  `h0 = 0.01*hu` failed (ConvergenceFailed) on perturbed Robertson at t = 4e5, where the warm BDF
  start succeeded. Large carried h0 values with a cold Adams start can fail in this way.
- The warm start does not always save nfe. On perturbed restarts the cold start was often cheaper
  (Jacobian cost). Whether radkin keeps it is for the radkin-side measurement, per the grill's
  keep rule.

## Commits

Filled by ticket-complete.js.

## Time Spent

Filled by ticket-complete.js.
