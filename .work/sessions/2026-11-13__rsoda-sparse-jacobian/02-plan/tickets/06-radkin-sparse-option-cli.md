# Ticket 06: radkin-sparse-option-cli

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `SolveOptions` gains `jacobian: JacobianMode { Dense, Sparse }` (default Dense); `use_jac = false` with Sparse selects `SparseFd`, otherwise `SparseUser`; the `OdeSystem` impl in `solver.rs` forwards the sparse callback to `jacobian_sparse`.
- [ ] The CLI accepts `--sparse-jacobian` and forwards it to the run (chunked and single-pulse drivers).
- [ ] The golden run tests pass in sparse mode (analytic and FD) at the existing rtol 1e-6, without loosening.
- [ ] `examples/bench.rs` reports dense vs sparse timings; README documents the option and flag.

**Files to Touch:**
- `src/solver.rs`
- `src/unit_cell.rs`
- `src/bin/radkin.rs`
- `tests/` (golden run tests parameterised over the Jacobian mode)
- `examples/bench.rs`
- `README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example bench
```

Expected:
all green; the golden runs pass for Dense and Sparse; bench prints dense and sparse timings.

**Notes:**

If a sparse golden run misses rtol 1e-6, compare its step sequence with the dense one first (`nje`, `nfe`, method switches); a missing pattern entry shows up as a Jacobian mismatch in ticket-05's test, not here.
