# Ticket 02: rsoda-warm-start

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] An opt-in warm start on restart (e.g. `restart_warm` or a `RestartMode` argument) keeps the previous solve's method (Adams/BDF) and order; the default `restart` is unchanged.
- [ ] If starting above order 1 without history proves invalid (error-test failures or loss of accuracy in a test), the warm start keeps only the method, and the doc comment says so.
- [ ] A test shows a warm-started solve stays within tolerance of a fresh solve on a stiff problem; rsoda checks are clean; the branch is pushed.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/src/lsoda.rs`
- `${DEV_DIR}/RUST/rsoda/tests/`

**Verification Step:**

Run:
```bash
cd "$DEV_DIR/RUST/rsoda" && cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features
```

Expected:
All checks clean, all tests pass.

**Notes:**

Same branch `feat/lsoda-restart`. Record in the commit message whether the order is carried or only the method, and why.
