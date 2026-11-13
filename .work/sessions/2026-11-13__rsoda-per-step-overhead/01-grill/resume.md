# Session: rsoda-per-step-overhead

**Date:** 2026-11-13T08:31:51.721Z
**Status:** Grill phase complete

## Problem Statement

We do not know how a radkin solve's time splits between the chemistry kernel (derive, Jacobian) and rsoda (LU factor, solve, bookkeeping), so we cannot target solver overhead before Phase 2 multiplies the per-RHS cost by the number of voxels. This session measures the time breakdown, then keeps the optimizations that the data justifies.

## Context & Constraints

- **Current behavior:** Labarbe run (1e-7..10 s, dense): solve 7.5 ms, nfev 4298, nst 2567, so a per-RHS cost of 1.75 us. With n_dyn 22, the kernel alone measures 0.57 us per `derive` and 2.0 us per dense `jacobian`, so the kernel is about 1/3 of the solve. A grep of rsoda's step loop (stoda, prja, solsy, linalg) finds no per-step allocations. radkin's `solve` stores every step (`y.to_vec()` per step) when t_eval is None.
- **Pain point:** The remaining ~2/3 (rsoda factor/solve/bookkeeping, output storage) is not broken down, so there is no data to choose optimizations from.
- **Dependencies:** rsoda is a git dependency pinned by rev: rsoda changes are made and pushed in the rsoda repo, then radkin bumps the `rev` (no committed `path` dependency). The golden tests of both repos (rsoda golden at 1e-12; radkin RHS/Jacobian/G-values at 1e-12, runs at rtol 1e-6) must pass at their current tolerances. CI (fmt, clippy -D warnings, tests) stays clean on Ubuntu and Windows in both repos.
- **Tech stack:** Rust 2024, rsoda (LSODA port, dense/banded/faer sparse linear solvers behind the LinearSolver trait), examples/bench.rs, tools/benchmark.py.

## Success Metrics

- The instrumented bench prints a time breakdown (kernel RHS, kernel Jacobian, factor, solve, remainder) for every case of the tools/benchmark.py matrix, dense and sparse, and docs/perf/ records the baseline and how to reproduce it.
- Every kept optimization cuts the median solve time by at least 5% on its target cases, with no case slower beyond noise. Rejected attempts are recorded in docs/perf/ with their numbers.
- All golden tests in rsoda and radkin pass at unchanged tolerances, and CI is green in both repos. radkin pins the new rsoda rev.

## Architecture & Approach

1. Instrument: rsoda gets an optional `timing` cargo feature (off by default, zero cost when off) that counts and times factor and solve in the linear solver. examples/bench.rs wraps the OdeSystem to time rhs/jac calls and prints the time breakdown for each case, dense and sparse.
2. Profile: run the tools/benchmark.py matrix through the instrumented bench and write docs/perf/ with the baseline breakdown and reproduction steps.
3. Optimize: one ticket per candidate that the data points to (e.g. kernel derive, dense LU/solve, rsoda bookkeeping, output storage). A change stays only if it meets the 5% rule. There is no global target: the session ends when the profile shows no cheap wins left.
4. rsoda changes are pushed to rsoda and radkin bumps the rev. docs/perf/ gets the final numbers, plus what was tried and dropped.

## Assumptions & Trade-offs

- Numerics: reordering floating-point operations (e.g. kernel summation order, LU loop order) is allowed as long as the golden tests pass at their current tolerances. The step heuristics (Jacobian reuse, step-size control) stay as in ODEPACK, to keep parity with scipy/radiopyo.
- Profiling uses counters and timers kept in the repos, not a sampling profiler (none is installed besides WPR), so the numbers repeat on CI and Linux.
- The cases are the existing benchmark matrix. A synthetic large-N (Phase-2-scale) system is out of scope.
- No radiopyo is needed for before/after numbers. Only radkin's own solve time is compared.
- The concrete optimization tickets depend on the baseline profile, so the plan drafts them from the measured breakdown.

## Open Questions

- Which optimizations are worth doing depends on the baseline breakdown measured in the first ticket.
- Timer overhead (Instant::now per call at ~0.5 us granularity) may distort the breakdown. The instrumented totals are checked against an uninstrumented run.

## Notes

- The idea is `.work/idea/rsoda-per-step-overhead.md`. Its "avoid per-step allocations" suggestion looks already done in rsoda (grep), and "check Jacobian reuse heuristics vs ODEPACK" is out of scope (heuristics stay as in ODEPACK).
- Glossary: added Time breakdown and Per-RHS cost.
