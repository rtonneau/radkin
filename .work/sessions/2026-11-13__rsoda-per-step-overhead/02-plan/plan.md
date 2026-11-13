# Implementation Plan

**Session:** rsoda-per-step-overhead
**Date:** 2026-11-13T08:40:20.085Z
**Estimated effort:** 1-2 days

## Strategy

Measure first, then optimize. Ticket 1 adds an optional `timing` feature to rsoda that counts and times the linear solver's factor and solve calls. Ticket 2 has radkin's bench use it, adds rhs/jac timers, prints the time breakdown for the benchmark matrix (dense and sparse) and writes the baseline to docs/perf/. Tickets 3-5 each go after one candidate the baseline points to: the kernel's derive/jacobian in radkin, the dense LU factor/solve in rsoda, and rsoda's per-step bookkeeping plus radkin's output storage. Each of them keeps a change only if the median solve time drops by at least 5% on its target cases with no case slower beyond noise. Otherwise the attempt is reverted and recorded in docs/perf/. Ticket 6 pins radkin to the final pushed rsoda rev and writes the final numbers.

rsoda work goes on a branch `perf/per-step-overhead` in the rsoda repo with its own PR. During development radkin uses a temporary local `[patch]` for rsoda, which is never committed. Committed radkin code pins a pushed rsoda rev.

## Tickets Overview

- **Ticket 1:** rsoda `timing` cargo feature: counters and timers for linear-solver factor and solve, exposed through a stats accessor, zero cost when off.
- **Ticket 2:** radkin bench prints the time breakdown (kernel RHS, kernel Jacobian, factor, solve, remainder) per case, dense and sparse. docs/perf/ records the baseline and how to reproduce it.
- **Ticket 3:** Optimize the radkin kernel (derive, dense and sparse jacobian) where the profile shows the cost, under the 5% rule.
- **Ticket 4:** Optimize rsoda's dense LU factor/solve (DGEFA/DGESL port) under the 5% rule.
- **Ticket 5:** Optimize rsoda's per-step bookkeeping (stoda/lsoda loop, norms, Nordsieck updates) and radkin's per-step output storage under the 5% rule.
- **Ticket 6:** Push rsoda, pin radkin to the new rev, rerun the matrix, and write the final breakdown plus what was tried and dropped in docs/perf/.

## Sequencing Rationale

The instrumentation (1) must exist before the profile (2), and the profile picks which optimizations are worth the effort (3-5). Tickets 3-5 touch separate code (radkin kernel, rsoda linalg, rsoda step loop and radkin output), so each is measured against the same baseline without interference. Ticket 6 depends on all rsoda changes being final.

## Risks & Mitigation

- **Risk:** Timer overhead (Instant::now around ~0.5 us calls) distorts the breakdown. → **Mitigation:** Report the instrumented total next to an uninstrumented run of the same case and give the overhead in docs/perf/. Optimization decisions use uninstrumented median solve times only.
- **Risk:** Run-to-run noise on Windows hides a 5% effect. → **Mitigation:** Use --repeat >= 10 with medians, run before/after back to back on the same machine, and repeat a borderline result before keeping it.
- **Risk:** Reordered floating-point operations change step counts and push a run past a golden tolerance. → **Mitigation:** Run both repos' golden tests after each change. If one fails, revert to the original order for that code. Never loosen a tolerance.
- **Risk:** radkin's committed Cargo.toml points at an unpushed rsoda rev. → **Mitigation:** Ticket 6 pushes rsoda first, then bumps the rev. The `[patch]` used during development is removed before each radkin commit.

## Assumptions

- `DEV_DIR` resolves to `D:\DEV`. rsoda is checked out at `${DEV_DIR}/RUST/rsoda`, on a clean `main`.
- The benchmark matrix in tools/benchmark.py can be driven through examples/bench.rs without radiopyo.
- Step heuristics stay as in ODEPACK. Only the cost per operation is reduced.
