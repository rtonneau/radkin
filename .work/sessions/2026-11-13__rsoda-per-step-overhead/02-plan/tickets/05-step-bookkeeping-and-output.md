# Ticket 05: step-bookkeeping-and-output

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Based on the breakdown's remainder, hot spots in rsoda's step loop (stoda corrector loop, error norms vmnorm, Nordsieck rescale and prediction, ewt updates) are optimized without changing step heuristics.
- [ ] radkin's per-step output storage in `solver::solve` (one `Vec` per stored step) is optimized only if measurable, without changing the public `Solution` layout unless the gain meets the 5% rule and every caller is updated.
- [ ] Each change is kept only under the 5% rule. Otherwise it is reverted.
- [ ] Golden tests in both repos pass unchanged, and docs/perf/README.md records each attempt.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/src/stoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/lsoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/norms.rs`
- `src/solver.rs`
- `docs/perf/README.md`

**Verification Step:**

Run:
```bash
cd $DEV_DIR/RUST/rsoda && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features && cd $DEV_DIR/RUST/radkin && cargo clippy --release --all-targets -- -D warnings && cargo test --release && python -I tools/breakdown.py --repeat 10
```

Expected:
All tests pass in both repos, and the breakdown shows a remainder drop for each kept change (or a recorded no-op).

**Notes:**

The step counts (nst, nfe, nje) must stay as before on every case, since the heuristics are unchanged. A changed count signals a numerics change that needs scrutiny against the golden tests.
