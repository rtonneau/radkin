# Ticket 01: rsoda-linear-solver-trait

**Status:** ✅ Done

## Local Test Result

Run in `D:/DEV/RUST/rsoda` (branch `feat/sparse-jacobian`):
`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → rc=0; 107 unit + 10 analytic + 1 golden + 1 doc tests pass; no tolerance edits.

Extra bit-identity checks:
- Robertson probe (DenseUser, DenseFd, BandedUser{2,2}, BandedFd{2,2}; 12 output times each): `y` bits, `nst`, `nfe`, `nje` identical between 312b742 and this branch (`diff` empty).
- radkin `cargo test --release` against the local rsoda via `--config patch...path` (not committed): all pass.

## Review Notes

Base: rsoda `312b74205f0d6061037f5a642522bb2a35ae6a27` (main), clean tree, new branch `feat/sparse-jacobian`. radkin base: session branch `feat/rsoda-sparse-jacobian` (only gps session files changed).

Acceptance Criteria:
- Public `LinearSolver` trait (`src/linsol.rs`): `prepare(kind, &mut dyn OdeSystem, &PrjaInput, y, ftem, nfe, nje) -> PrjaOutput` (DPRJA) and `solve(h, el0, x)` (DSOLSY). Supertrait `Send`, so `Lsoda` stays `Send`. Re-exported with `MatrixSolver`, `PrjaInput`, `PrjaOutput`, `SingularMatrix`.
- `MatrixSolver` wraps `prja` and `IterationMatrix::solve` verbatim (dense, banded; diagonal DSOLSY code untouched), so the arithmetic is unchanged.
- `State.wm: Box<dyn LinearSolver>` plus `wm_kind`. stoda calls `st.wm.prepare` / `st.wm.solve`. `lsoda::match_matrix_storage` reallocates via `same_storage`/`take_solver`/`release_solver`; `start()` recycles a custom solver. The kind-to-solver mapping lives in `linsol::builtin_solver`/`same_storage`.
- Custom solvers: `JacobianKind::Custom` + `Lsoda::with_linear_solver(system, options, t0, y0, Box<dyn LinearSolver>)`; `InputError::MissingLinearSolver` for `new` or `set_jacobian(Custom)` without one. Tests: a counting wrapper around the dense solver is bit-identical to DenseUser on Robertson (same y, nst, nfe, nje; prepare called nje times); kept across `set_jacobian`; missing-solver errors.
- Existing tests unchanged except internal asserts on `wm.storage` → `wm_kind` (same meaning).

Design notes: `MatrixStorage::for_jacobian` now returns `Option` (`None` for Custom). An unsized system is passed to the solver through a `DynSystem` adapter (crate-private). `State` no longer derives `Clone`/`PartialEq` (they were unused).

Extra commit `a9d58f1` (formatting only): rsoda at 312b742 was not `cargo fmt`-clean under any config. `rustfmt.toml` (`max_width = 120`, `use_small_heuristics = "Max"`) gave the smallest diff from the existing style; the whole crate was formatted in that separate commit so the ticket diff holds only real changes. Pre-existing clippy errors (Fortran-style index loops, NaN-aware `!(x >= 0)`, LINPACK `>= y + 1` checks, argument counts) are allowed crate-wide in `[lints.clippy]` in Cargo.toml instead of rewriting the ported code.

Review: inline (no independent reviewer spawned; the change was checked by the bit-identity probe plus the full rsoda and radkin suites).

## Blockers / Challenges

- rsoda had no fmt/clippy-clean baseline (see Review Notes); fixed with a formatting-only commit and crate-level lint allows.

## Commits

- c9d0010 chore(rsoda): ticket 01 LinearSolver trait at rsoda f6b24eb (feat/sparse-jacobian)

## Time Spent

11m (ticket-start.js to ticket-complete.js)
