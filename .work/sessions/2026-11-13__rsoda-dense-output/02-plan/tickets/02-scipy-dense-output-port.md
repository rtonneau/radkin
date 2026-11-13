# Ticket 02: scipy-dense-output-port

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] A pure function (e.g. `src/dense.rs`) ports scipy `_dense_output_impl` + `LsodaDenseOutput._call_impl`: columns `0..=nqu`, last column times `(hcur/hu)^nqu` when `nqcur < nqu`, `x = (t - t_step)/hcur`, `y_i = sum_j yh[i, j] * x.powf(j)`; no allocation in the hot path (caller buffers).
- [ ] `tools/export_golden.py` exports scipy `LsodaDenseOutput` reference values for synthetic `yh` cases (order kept, order decreased, negative x, x = 0) into `tests/golden/`.
- [ ] A Rust golden test compares at 1e-12 relative.

**Files to Touch:**
- `src/dense.rs` (new) and `src/lib.rs`
- `tools/export_golden.py`
- `tests/golden/dense_output.json` (new)
- `tests/` golden test file

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks pass, including the new dense-output golden test.

**Notes:**

numpy `float ** int array` is C `pow`, so use `powf`, not `powi`. Mirror the rescale exactly: `yh[:, -1] *= (h / rwork[10]) ** order` with `h = HCUR`, `rwork[10] = HU`.
