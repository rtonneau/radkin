# Ticket 02: rsoda-sparsity-pattern-coloring

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] A public CSC `SparsityPattern` type (n, col_ptr, row_idx) with validation (sorted, in range, no duplicates) returning `InputError` variants; it always includes the diagonal (added if missing).
- [ ] Greedy Curtis-Powell-Reid column coloring: columns in the same group share no nonzero row; tests check validity on random patterns and that a tridiagonal pattern gets 3 colors.
- [ ] Colored-FD Jacobian assembly into CSC values using DPRJA's increment formula per column (`r = max(srur*|y_j|, r0/ewt_j)`), one `rhs` call per color; on a dense-pattern problem it reproduces the dense FD Jacobian to roundoff.
- [ ] A weighted norm over CSC values with the same definition as `fnorm` (equal to `fnorm` of the densified matrix in tests).
- [ ] No new dependencies; all existing tests still pass.

**Files to Touch:**
- `D:/DEV/RUST/rsoda/src/sparse.rs` (new; pattern, coloring, assembly)
- `D:/DEV/RUST/rsoda/src/norms.rs`
- `D:/DEV/RUST/rsoda/src/error.rs`
- `D:/DEV/RUST/rsoda/src/lib.rs`

**Verification Step:**

Run:
```bash
cd D:/DEV/RUST/rsoda && cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
all green; new unit tests for pattern validation, coloring validity, colored FD vs dense FD, and the sparse norm pass.

**Notes:**

Keep these helpers usable without the `sparse` feature: they are pure and FD-only users or custom solvers may want them. Match DPRJA's `fac`/`r0` computation (`r0 = 1000*|h|*UROUND*n*fac`) exactly as the dense path does.
