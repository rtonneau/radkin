# Ticket 04: rsoda-dense-lu

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Based on the breakdown, rsoda's dense factor and solve (DGEFA/DGESL port and iteration-matrix assembly in linalg.rs/linsol.rs/prja.rs) are optimized where they matter (e.g. column-major inner loops that vectorize, fused assembly and factor, bounds-check removal).
- [ ] A change is kept only under the 5% rule (median uninstrumented solve time of the radkin bench matrix, dense). Otherwise it is reverted.
- [ ] rsoda's golden (1e-12), analytic and sparse tests pass unchanged, and radkin's golden tests pass with the patched rsoda.
- [ ] docs/perf/README.md records each attempt with numbers.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/src/linalg.rs`
- `${DEV_DIR}/RUST/rsoda/src/linsol.rs`
- `${DEV_DIR}/RUST/rsoda/src/prja.rs`
- `docs/perf/README.md`

**Verification Step:**

Run:
```bash
cd $DEV_DIR/RUST/rsoda && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release --all-features && cd $DEV_DIR/RUST/radkin && cargo test --release && python -I tools/breakdown.py --repeat 10
```

Expected:
All tests pass in both repos, and the breakdown shows a factor/solve share drop for each kept change (or a recorded no-op).

**Notes:**

Keep LINPACK pivoting semantics (partial pivoting, same pivot choice), since a different pivot changes the results beyond reordering. Banded LU is out of scope unless trivially shared. radkin is tested through the uncommitted local `[patch]`.
