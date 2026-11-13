# Ticket 03: rsoda-faer-sparse-solver

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] Optional cargo feature `sparse` adds a `faer` dependency; with the feature off, `cargo tree` shows no dependencies.
- [ ] A `SparseLu` `LinearSolver`: symbolic analysis once per pattern, numeric refactorization of `P = I - h*el0*J` at each Jacobian update, a singular `P` reported through `PrjaOutput::singular`, and `pdnorm` from the ticket-02 sparse norm.
- [ ] `JacobianKind::SparseUser { pattern }` and `JacobianKind::SparseFd { pattern }` (behind the feature); `OdeSystem` gains a sparse Jacobian callback writing CSC values in pattern order (default `MissingJacobian`), and `ClosureSystem` supports it; option validation covers the new kinds.
- [ ] A Brusselator 2D reaction-diffusion test (a grid of at least 16x16 cells, 2 species, BDF-dominated): SparseUser and SparseFd match BandedUser at rtol 1e-8 to within 1e-6 relative on the final state; SparseFd uses `ncolors` RHS calls per Jacobian (checked via `nfe`/`nje`).
- [ ] Method switching (Adams to BDF and back) works with the sparse kinds (a test where `mused` changes).

**Files to Touch:**
- `D:/DEV/RUST/rsoda/Cargo.toml`
- `D:/DEV/RUST/rsoda/src/sparse.rs`
- `D:/DEV/RUST/rsoda/src/options.rs`
- `D:/DEV/RUST/rsoda/src/system.rs`
- `D:/DEV/RUST/rsoda/src/lsoda.rs`
- `D:/DEV/RUST/rsoda/tests/sparse.rs` (new)

**Verification Step:**

Run:
```bash
cd D:/DEV/RUST/rsoda && cargo clippy --release --all-targets --features sparse -- -D warnings && cargo test --release --features sparse && cargo test --release && cargo tree -e normal --depth 1
```

Expected:
both test runs pass; the Brusselator tests pass; `cargo tree` (feature off) lists only `rsoda`.

**Notes:**

Read faer's sparse LU source for the chosen version before coding (symbolic/numeric split, partial pivoting, the error on a singular matrix); confirm, don't guess. Keep all faer usage in one module. Assemble `P`'s values into a buffer preallocated from the pattern (no per-step allocation in the assembly; record whether faer's numeric factorization allocates). If the 1e-6 Brusselator bar fails, find out why (pattern, norm, scaling) instead of loosening it.
