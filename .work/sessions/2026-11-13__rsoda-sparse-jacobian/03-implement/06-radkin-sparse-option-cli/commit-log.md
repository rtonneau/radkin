# Ticket 06: radkin-sparse-option-cli

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example bench -- data/config_Labarbe.toml --t-span 1e-7 10 --repeat 3 --jacobian both` → rc=0, 14 test binaries ok.

- `tests/golden_runs.rs`: all four golden tests (labarbe_run, labarbe_curve, chunked_pulsed_run, paper_points) pass for dense, sparse and sparse-FD at the existing tolerances (rtol 1e-6, unchanged).
- `tests/cli.rs::sparse_jacobian_flag_matches_dense_run`: the final row of `--sparse-jacobian` matches the dense run within 1e-6.
- Bench, Labarbe run to 10 s: dense solve ≈ 7.3 to 7.8 ms, sparse ≈ 7.5 to 8.6 ms (nst 2567 vs 2573). Labarbe chunked to 1 s: dense ≈ 6.8 ms, sparse ≈ 7.5 ms. At 15 to 22 species, sparse is about the same speed as dense, as the design expected.

## Review Notes

Base: radkin `ae23dfc`.

- `solver::JacobianMode { Dense (default), Sparse }` and `SolveOptions.jacobian`. Kind mapping: (Dense, jac) → DenseUser, (Dense, no jac) → DenseFd, (Sparse, jac) → SparseUser, (Sparse, no jac) → SparseFd. `ChemSystem` holds `Option<SparseJacobian>` (built in `solve` only in sparse mode) and implements `sparsity_pattern` / `jac_sparse` / `has_jac_sparse`.
- `unit_cell::RunOptions.jacobian` (default Dense) is passed to `SolveOptions` for run, chunked and single-pulse (`ChunkOptions.run`).
- CLI `--sparse-jacobian` sets it. README documents the flag and the bench option.
- `examples/bench.rs`: `--jacobian dense|sparse|both`. With the default, the JSON keys are unchanged for `tools/benchmark.py`; `both` adds a `"sparse"` object.
- Note: the `SparseJacobian` (pattern, coloring, COLAMD) is rebuilt for each `solve` call, so once per chunk in chunked runs. That's negligible at this size; if it shows up in Phase 2, cache it next to `SimMatrices`.

Review: inline.

## Blockers / Challenges

None.

## Commits

- fbd0f3d feat(solver): JacobianMode::Sparse option, --sparse-jacobian flag, sparse golden runs and bench (ticket 06)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
