# Ticket 01: rsoda-timing-feature

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] rsoda has a `timing` cargo feature, off by default. With it on, the solver counts and times every linear-solver factor (assemble + factor) and solve call. Without it, there are no timers in the code path (cfg-gated) and results are bit-identical.
- [ ] A public accessor (e.g. `Lsoda::timing() -> TimingStats { n_factor, t_factor, n_solve, t_solve }`) exists only under the feature and is documented.
- [ ] A test under the feature checks that the counts match the solver's own counters (factor count equals the number of iteration-matrix evaluations; solves are > 0).
- [ ] rsoda CI runs clippy and tests with the feature on as well (added to the existing feature matrix).
- [ ] Work is on rsoda branch `perf/per-step-overhead`, committed.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/Cargo.toml`
- `${DEV_DIR}/RUST/rsoda/src/lsoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/stoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/lib.rs`
- `${DEV_DIR}/RUST/rsoda/.github/workflows/ci.yml`
- `${DEV_DIR}/RUST/rsoda/README.md`

**Verification Step:**

Run:
```bash
cd $DEV_DIR/RUST/rsoda && cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release && cargo test --release --features timing
```

Expected:
fmt and clippy clean; all tests pass with and without the feature.

**Notes:**

Find where stoda calls the LinearSolver's assemble/factor and solve (prja/solsy), and wrap only those call sites, so every solver implementation (dense, banded, sparse, custom) is covered. Use std::time::Instant. Do not change any numerics.
