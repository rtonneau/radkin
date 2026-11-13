# Ticket 03: rsoda-faer-sparse-solver

**Status:** ✅ Done

## Local Test Result

In `D:/DEV/RUST/rsoda`: `cargo fmt --check && cargo clippy --release --all-targets --features sparse -- -D warnings && cargo clippy --release --all-targets -- -D warnings && cargo test --release --features sparse && cargo test --release && cargo tree -e normal --depth 1` → rc=0.
- With `sparse`: 118 unit, 10 analytic, 1 golden, 4 sparse, 1 doc tests. Without it: 117 unit (including the feature-disabled check), 10, 1, 0 sparse (cfg-gated), 1.
- `cargo tree` with the feature off: `rsoda` only.
- Brusselator 2D, m = 16 (n = 512), rtol 1e-8, atol 1e-10, t = 10: SparseUser vs BandedUser{32,32} max rel. diff 4.6e-15; SparseFd 1.2e-15; same nst (791) and nje (36). Bar: 1e-6.
- Dense/band bit-identity probe vs 312b742 re-run with the feature on: identical.

## Review Notes

Base: rsoda `baa7fe5`; commit `81638e2`.

Acceptance Criteria:
- `sparse = ["dep:faer"]` with faer 0.24 (`default-features = false`, `std` and `sparse-linalg`). Feature off → no dependencies.
- `SparseLu` (`src/sparse_lu.rs`): COLAMD ordering once in `new`; each `prepare` assembles `-hl0*J` (analytic `jac_sparse` + `scale_jacobian`, or `fd_jacobian`), computes `pdnorm` with `finish_iteration_matrix`, adds I, and refactors with `faer::sparse::linalg::lu::simplicial::factorize_simplicial_numeric_lu` (partial pivoting by largest magnitude, factor vectors and `MemBuffer` reused). Singular: faer's `SymbolicSingular` or a zero U pivot. I read faer's source: a numerically zero pivot is not an error there, since it takes `recip()` of 0, so rsoda checks the U diagonal itself.
- `JacobianKind::SparseUser` / `SparseFd`, with the pattern from `OdeSystem::sparsity_pattern` (new, default `None`), `jac_sparse`, `has_jac_sparse`; `ClosureSystem::with_sparsity` / `with_sparse_jac`. Validation in `OptionsBuilder::build`: `MissingSparsityPattern`, `SparsityPatternSize`, `MissingJacobian`, `SparseFeatureDisabled`.
- Brusselator tests: see above. FD calls: a counting custom solver wrapping `SparseLu` sees exactly `ncolors` f calls per Jacobian (8 colors for n = 128); the same run via `SparseFd` is bit-identical.
- Method switching: Robertson with SparseUser and SparseFd switches Adams → BDF (`on_method_switch`) and stays within 1e-6 of DenseUser at six output times.

Deviations from the ticket text (implementation choices):
- The kinds carry no `pattern` field; the system provides it. This keeps `JacobianKind` `Copy` and the pattern next to the Jacobian callback that fills it.
- The variants are always present rather than cfg-gated. Without the feature they fail validation with `InputError::SparseFeatureDisabled`, which avoids feature-dependent enum shapes (downstream exhaustive matches would break under feature unification).
- Only faer's simplicial LU is used (not its supernodal variant). Supernodal would only matter for large fill; the trait makes adding it later local.

Review: inline.

## Blockers / Challenges

- At the first diffusion coefficient (0.002/dx²) the Brusselator was not stiff, so LSODA stayed in Adams. The test now uses 0.1/dx² to t = 10, which is BDF-dominated.

## Commits

- 90a301c chore(rsoda): ticket 03 SparseLu and sparse Jacobian kinds at rsoda 81638e2 (feat/sparse-jacobian)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
