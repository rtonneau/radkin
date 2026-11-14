# Ticket 01: rsoda-restart

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `Lsoda::restart(&mut self, t0: f64, y0: &[f64], h0: f64) -> Result<(), LsodaError>` exists in rsoda and keeps the system, the `State` buffers and the linear solver (no reallocation of the n x n or sparse storage).
- [ ] After `restart`, integrating gives bitwise the same y, nst, nfe and nje as a fresh `Lsoda::new` with the same options, on several states, in both the dense and sparse Jacobian modes.
- [ ] rsoda's fmt, clippy (`-D warnings`) and tests are clean; the commit is pushed on the rsoda branch `feat/lsoda-restart`.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/src/lsoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/state.rs`
- `${DEV_DIR}/RUST/rsoda/tests/` (a new restart test)

**Verification Step:**

Run:
```bash
cd "$DEV_DIR/RUST/rsoda" && cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features restart
```

Expected:
fmt and clippy are silent; the restart tests pass.

**Notes:**

Work in the rsoda repo, on a new branch `feat/lsoda-restart` from the commit radkin pins (`f5d1312`), not on the checked-out `perf/per-step-overhead`. Read what `Lsoda::new`/`with_linear_solver` and `State::with_solver` initialise, and reset every one of those fields. The rebuild path around `lsoda.rs:607` is a hint for what has to be reset. Push the branch to origin (not to `main`).
