# Ticket 05: radkin-sparse-kernel-jacobian

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `Cargo.toml` pins rsoda to the ticket-04 head SHA with `features = ["sparse"]` (git `rev`, no `path`).
- [ ] `SimMatrices` (or a new kernel item) builds the CSC sparsity pattern of the `n_dyn x n_dyn` Jacobian once: the union over K reactions (reactant columns x stoichiometry rows, plus switch-derivative columns from `negative_index`), enzymatic reactions (enzyme/substrate/negative columns), the overshoot rows (row j, the columns where `grad` can be nonzero), and the diagonal.
- [ ] `jacobian_sparse(y, mat, ws, values)` writes CSC values in pattern order without allocating (uses `Workspace` buffers), mirroring `jacobian`'s arithmetic.
- [ ] A test densifies the sparse fill and matches `jacobian` at 1e-12 for every golden config, including states with negative components (overshoot paths).
- [ ] The existing radkin tests still pass on the new rsoda rev (dense path unchanged).

**Files to Touch:**
- `Cargo.toml`
- `Cargo.lock`
- `src/kernel.rs`
- `src/matrices.rs`
- `tests/` (new sparse Jacobian test file)

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
all green, including the new sparse vs dense Jacobian tests and the unchanged golden tests.

**Notes:**

The scatter in `jacobian` adds `deriv` entries into `pd[row + col*n]`; for CSC, precompute per reaction the value-slot indices so the fill is an index lookup (no search) in the hot loop. Mind that `deriv` may contain repeated columns (accumulate). Column/row indices are for dynamic species only (`n_dyn`), as in the dense layout.
