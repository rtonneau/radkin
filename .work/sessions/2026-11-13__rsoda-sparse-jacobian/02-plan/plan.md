# Implementation Plan

**Session:** rsoda-sparse-jacobian
**Date:** 2026-11-13T07:37:05.435Z
**Estimated effort:** 3-4 days

## Strategy

The work spans two repositories. Tickets 01-04 are implemented in the rsoda checkout (`${DEV_DIR}/RUST/rsoda`) on its own branch `feat/sparse-jacobian`, pushed with a PR opened there. Tickets 05-06 are implemented in radkin on this session's branch, with `Cargo.toml` pinned to the pushed rsoda commit (a git `rev`, never a committed `path`; a temporary `[patch]` is fine while iterating). The order starts with a refactor that changes no behaviour and is proven bit-identical by rsoda's golden tests, then adds pure sparse building blocks (pattern, coloring, assembly, norm), then the faer-backed solver and the user-facing Jacobian kinds, then CI/docs. radkin consumes it last.

## Tickets Overview

- **Ticket 1:** rsoda: public `LinearSolver` trait; dense and banded become implementations; `prja`/`solsy`/`lsoda`/`stoda` go through it; bit-identical results.
- **Ticket 2:** rsoda: CSC `SparsityPattern`, greedy Curtis-Powell-Reid column coloring, colored-FD and analytic sparse Jacobian assembly, and the weighted sparse norm (zero dependencies).
- **Ticket 3:** rsoda: optional `sparse` feature with a faer sparse LU `LinearSolver` (symbolic analysis once, numeric refactorization per update), `JacobianKind::{SparseUser, SparseFd}`, a sparse `OdeSystem` Jacobian callback, and Brusselator 2D tests.
- **Ticket 4:** rsoda: CI (fmt/clippy/test on Ubuntu and Windows, feature off and on), README/docs, version bump; push the branch and open the rsoda PR.
- **Ticket 5:** radkin: bump rsoda `rev` with the `sparse` feature; kernel sparsity pattern and non-allocating sparse Jacobian fill, checked against the dense Jacobian at 1e-12.
- **Ticket 6:** radkin: `SolveOptions.jacobian: Dense | Sparse`, `--sparse-jacobian` CLI flag, golden runs in sparse mode at rtol 1e-6, benchmark entry, README.

## Sequencing Rationale

Ticket 1 isolates the risky refactor so the golden tests can prove it changes nothing before any new numerics exist. Ticket 2's pure functions are unit-testable without faer and are needed by both the feature-gated solver and the FD path. Ticket 3 depends on 1 (trait) and 2 (pattern/coloring/norm). Ticket 4 needs the full API to document and to give CI something to cover, and it publishes the commit radkin pins. radkin's Ticket 5 needs the pushed rsoda API. Ticket 6 needs the kernel fill from Ticket 5.

## Risks & Mitigation

- **Risk:** the trait refactor subtly changes dense/banded arithmetic (operation order, scaling) → **Mitigation:** move the existing code verbatim into the implementations; rsoda golden tests and radkin golden tests must pass unchanged before Ticket 2 starts.
- **Risk:** the faer sparse LU API differs from expectations (version churn, pivoting control, symbolic reuse) → **Mitigation:** read faer's source/docs for the pinned version first; wrap it behind the trait so only one file depends on it.
- **Risk:** sparse pivoting changes the LSODA step sequence, so sparse runs differ from dense by more than roundoff → **Mitigation:** the acceptance bar is radkin's existing golden rtol 1e-6, which must not be loosened; on a mismatch, investigate (pattern missing an entry, wrong norm) rather than relax.
- **Risk:** the kernel sparsity pattern misses entries (overshoot terms, enzymatic, switch derivatives), so values written outside the pattern are dropped → **Mitigation:** derive the pattern as a superset from the index tables; test that densifying the sparse fill equals the dense Jacobian at 1e-12 on every golden config, including the negative-concentration (overshoot) paths.
- **Risk:** rsoda's PR merge produces a different commit from the one radkin pins → **Mitigation:** pin radkin to the pushed rsoda branch head; after the rsoda PR merges, a follow-up bumps `rev` to the merge commit (noted in finish).

## Assumptions

- `DEV_DIR` = `D:\DEV`; the rsoda checkout at `D:\DEV\RUST\rsoda` is clean on `main` at `312b742`, and pushing a branch and opening a PR on rtonneau/rsoda is covered by the user's `/gps auto` request.
- rsoda has no CI workflow yet; Ticket 4 adds one.
- Breaking rsoda API changes are acceptable (0.x, git dependency only); `JacobianKind` gains variants and `IterationMatrix`/`MatrixStorage` may be replaced by the trait.
- The iteration matrix always needs the diagonal, so the sparse pattern is the union of the user pattern and the diagonal.
- Dense stays radkin's default; sparse is opt-in and is not expected to be faster for the 20-40 species single cell.
