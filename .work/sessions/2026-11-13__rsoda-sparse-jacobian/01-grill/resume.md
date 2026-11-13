# Session: rsoda-sparse-jacobian

**Date:** 2026-11-13T07:30:24.587Z
**Status:** Grill phase complete

## Problem Statement

rsoda only factors the iteration matrix `P = I - h*el0*J` as dense or banded LU. Large stiff systems (for example a Phase-2 voxel grid with n_voxels x n_species unknowns) are intractable that way. Give rsoda a pluggable linear-solver layer with a sparse LU path, and let radkin opt into a sparse kernel Jacobian, without committing Phase 2 to operator splitting or a monolithic solve.

## Context & Constraints

- **Current behavior:** rsoda `MatrixStorage::{Dense, Diagonal, Banded}` in `solsy.rs`; `prja.rs` forms J (user or FD), scales it, computes `pdnorm` with `fnorm`/`bnorm`, adds the identity and factors with `dgefa`/`dgbfa`. `JacobianKind::{DenseUser, DenseFd, BandedUser, BandedFd}`. `OdeSystem::jac` writes a dense or band column-major `pd`. radkin `src/solver.rs` always uses `DenseUser` (or `DenseFd`) and `src/kernel.rs::jacobian` fills a dense `n_dyn x n_dyn` matrix.
- **Pain point:** dense/banded LU cost and memory grow as n^3 / n*bw^2, so a voxelised monolithic solve is impossible. A user-supplied solver cannot be plugged in.
- **Dependencies:** rsoda dense/banded results must stay bit-identical (rsoda golden tests). rsoda keeps zero dependencies when the new feature is off. radkin pins rsoda by git `rev`: the change is pushed in rsoda first, then radkin bumps `rev` (no committed `path` dependency). radkin golden tests (RHS/Jacobian/G-values at 1e-12, runs at rtol 1e-6) must not be loosened.
- **Tech stack:** Rust 2024; rsoda (LSODA port, `src/prja.rs`, `src/solsy.rs`, `src/norms.rs`, `src/options.rs`, `src/system.rs`, `src/stoda.rs`, `src/lsoda.rs`); faer sparse LU behind an optional rsoda `sparse` cargo feature; radkin `src/kernel.rs`, `src/matrices.rs`, `src/solver.rs`, `src/bin/radkin.rs`, `examples/bench.rs`.

## Success Metrics

- rsoda: all existing tests and golden data pass unchanged (bit-identical dense/banded) with the `sparse` feature both off and on; CI covers both.
- rsoda: on a 2D Brusselator reaction-diffusion grid, the sparse path (analytic and colored-FD) matches the banded solution within a tight tolerance, and colored FD uses one RHS call per color per Jacobian.
- radkin: with `SolveOptions.jacobian = Sparse` (and `--sparse-jacobian`), the golden run tests pass at rtol 1e-6, the sparse kernel Jacobian matches the dense one at 1e-12, and the kernel does not allocate; `examples/bench.rs` reports sparse vs dense timings. fmt/clippy/test stay clean on both repos.

## Architecture & Approach

rsoda:
1. A public `LinearSolver` trait that owns the iteration-matrix storage: assemble `P` from J (scaled by `-h*el0`), report the weighted norm of J (`pdnorm`), factor, and solve. Dense and banded become implementations that keep the exact current arithmetic. It is a public extension point, so users can plug their own solvers (e.g. Krylov, block-diagonal) later.
2. Sparse Jacobian kinds: the user gives a fixed CSC sparsity pattern. The analytic variant fills values through a new sparse Jacobian callback on `OdeSystem`. The FD variant uses greedy Curtis-Powell-Reid column coloring computed by rsoda from the pattern.
3. An optional `sparse` cargo feature with a faer-based sparse LU: one symbolic factorization for the fixed pattern, then a numeric refactorization at each Jacobian update. The weighted norm is computed from the CSC values with the same definition as `fnorm`.
4. Tests: Brusselator 2D grid comparing sparse with banded and sparse-FD with sparse-analytic; CI matrix with the feature off and on. A version bump with breaking API changes is acceptable (0.x, git dependency).

radkin:
5. Bump the rsoda `rev` and enable the `sparse` feature.
6. Kernel: derive the sparsity pattern from reactant/stoichiometry/enzyme indices plus the overshoot terms, and a non-allocating sparse Jacobian fill writing CSC values.
7. `SolveOptions.jacobian: Dense | Sparse` (Dense by default) and a `--sparse-jacobian` CLI flag. TOML configs are unchanged.
8. Golden runs in sparse mode, a sparse vs dense Jacobian check, and a benchmark entry.

## Assumptions & Trade-offs

- Phase 2's choice between operator splitting and a monolithic solve stays open. This is a general rsoda capability (ADR 0001).
- Not doing: a Krylov/LSODKR port, an automatic dense/sparse choice by size, a voxel prototype, or TOML options.
- The sparse LU's pivoting differs from dense, so sparse results agree with dense to solver tolerance, not bitwise. The golden rtol 1e-6 is the acceptance bar.
- For the ~20-40 species single cell, sparse may be no faster than dense. Dense stays the default.
- faer is a sizeable dependency, but it is opt-in through the feature.

## Open Questions

- The exact faer version and sparse LU API (simplicial vs supernodal, pivoting control) are settled during implementation.
- The precise trait signature (how FD Jacobian evaluation and `OdeSystem` interact with the solver) is settled in the plan.

## Notes

- Glossary: Iteration matrix, Linear solver, Sparsity pattern, Column coloring (`.work/GLOSSARY.md`). ADR 0001: pluggable linear solver with an optional faer sparse feature.
- Seeded from the Idea `.work/idea/rsoda-sparse-jacobian.md`. Related Ideas: phase2-voxel-diffusion, rsoda-per-step-overhead.
- The user approved the design and asked for `/gps auto`. The work is high-risk, so the drafted tickets still need explicit approval before `write:plan`.
