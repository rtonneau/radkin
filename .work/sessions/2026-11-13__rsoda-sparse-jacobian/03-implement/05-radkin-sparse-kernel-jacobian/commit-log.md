# Ticket 05: radkin-sparse-kernel-jacobian

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → rc=0; all 14 test binaries pass, including the unchanged golden tests (`tests/kernel.rs`, `tests/golden_runs.rs`) on rsoda 0.2.0.

New `tests/sparse_jacobian.rs`: the densified `jacobian_sparse` equals `jacobian` **exactly**, bit for bit (the spec asks for 1e-12):
- `config_Labarbe` (n = 22, nnz = 185, 1 enzymatic reaction), `advanced_water_radiolysis` and `CONFIG_PAPER_01` (n = 15, nnz = 164): 60 pseudo-random states each, over 9 decades with 0 %, 30 % and 100 % negative components (overshoot rows).
- The golden `synthetic` and `overshoot_state` states of labarbe and pure_water.

## Review Notes

Base: radkin `33b4d03` (session branch `feat/rsoda-sparse-jacobian`).

- `Cargo.toml`: `rsoda = { git, rev = "2f47645", features = ["sparse"] }` (PR rtonneau/rsoda#3 head). `Cargo.lock` updated (faer and its dependencies).
- `kernel::SparseJacobian::new(&SimMatrices)`: pattern = per K/enzymatic reaction its stoichiometry rows × dynamic columns (reactants or enzyme/substrate ∪ negative products), plus the overshoot row j × (slots ∪ slots_dj ∪ negative_without_j, or enzyme/substrate ∪ negative_without_j), plus the diagonal. It precomputes `Slots` (candidate columns and the value positions per column × row).
- `kernel::jacobian_sparse(y, mat, sj, ws, values)`: the same statements as `jacobian`, with `scatter_sparse` (a linear scan of the reaction's few candidate columns, then indexed writes) and the overshoot row through precomputed positions. No allocation (the `Workspace` buffers are reused).
- The candidate columns are a superset of every column the derivative lists can produce (mass-action slots, switch negatives, enzyme/substrate, overshoot gradient), so `expect` in `scatter_sparse` can't trigger. The exact-equality test across the overshoot paths confirms it.

Review: inline.

## Blockers / Challenges

None.

## Commits

- 08aa316 feat(kernel): sparse Jacobian pattern and fill on rsoda 0.2.0 (sparse) (ticket 05)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
