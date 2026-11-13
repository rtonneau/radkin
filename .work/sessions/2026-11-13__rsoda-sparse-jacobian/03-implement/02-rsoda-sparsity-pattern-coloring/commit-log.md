# Ticket 02: rsoda-sparsity-pattern-coloring

**Status:** ✅ Done

## Local Test Result

In `D:/DEV/RUST/rsoda`: `cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → rc=0 (116 unit, 10 analytic, 1 golden, 1 doc). `cargo tree -e normal --depth 1` → only `rsoda`.

## Review Notes

Base: rsoda `f6b24eb` (ticket 01) on `feat/sparse-jacobian`, clean tree. rsoda commit `baa7fe5`; radkin marker commit records it.

Acceptance Criteria:
- `SparsityPattern` (`src/sparse.rs`): `new` (CSC, validated, `InputError::InvalidSparsityPattern { reason }`), `from_entries`, `dense`, `banded`, accessors, `find`, `diag`. The diagonal is always inserted. Tests cover every validation branch and diagonal insertion.
- `ColumnColoring`: greedy CPR in column order via a row-wise transpose. Validity checked on 50 random patterns; tridiagonal → 3 colors, (ml=2, mu=1) → 4, dense → n singleton groups, empty → 1.
- `fd_jacobian`: DPRJA increments (`r0 = 1000|h|·UROUND·n·‖savf‖`, `r = max(srur|y_j|, r0/ewt_j)`), one `rhs` per color, `y` restored from `yh1` (also on `rhs` error, `nfe` then unchanged). On the dense pattern, `fd_jacobian` + `finish_iteration_matrix` + `dgefa` gives the same `pdnorm`, LU factors and pivots as `prja(DenseFd)`, bit for bit. Colored vs dense FD on a tridiagonal nonlinear system agree to 1e-9 relative, with zeros outside the pattern.
- `norms::csc_fnorm` equals `fnorm` of the densified matrix exactly (same summation order).
- `scale_jacobian` / `finish_iteration_matrix` (PDNORM, then identity on the diagonal) are ready for the analytic path of ticket 03.
- No dependencies added.

Review: inline.

## Blockers / Challenges

None.

## Commits

- 54b5605 chore(rsoda): ticket 02 sparsity pattern and coloring at rsoda baa7fe5 (feat/sparse-jacobian)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
