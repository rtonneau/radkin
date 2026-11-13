# Ticket 02: bench-time-breakdown-baseline

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo clippy --release --all-targets --features timing -- -D warnings && cargo test --release` → exit 0, 54 tests passed (log `02-verify.log`).
`python -I tools/breakdown.py --repeat 10 --out .scratch/.../02-baseline.json` → a table for all 9 cases x dense/sparse (18 rows), copied into docs/perf/README.md.

## Review Notes

Base: `52c0e472361e30571742fb43862074e06ad2520b`, clean apart from the gps session record.
- radkin feature `timing = ["rsoda/timing"]`, off by default. `solver::timing` (cfg-gated) is a thread-local `Breakdown` that rhs/jac/jac_sparse add to, plus rsoda's TimingStats at the end of each `solve`. `take()` returns and resets it, so chunked runs are summed. The default build has no timers, and the kernel still does not allocate.
- bench `--breakdown` errors without the feature, and adds per-repeat breakdowns to the dense and sparse lines.
- tools/breakdown.py builds a plain binary (target/) and a timing binary (target/timing), runs the benchmark.py CASES, and prints shares: rhs, jac, factor = prepare - jac, lin, rest, and the timer overhead. It reuses benchmark.py's helpers, which import pandas via compare_radiopyo. That works in the env used, but radiopyo itself is not needed.
- CI: added a clippy job with `--features timing`. AGENTS.md lists the tool.
- Baseline: rest 32-44% (largest everywhere), rhs 24-40%, lin 9-17%, factor 4-8% dense / 10-20% sparse, jac 4-8%. Timer overhead x1.02-1.13.

## Blockers / Challenges

A Python text-mode edit converted Cargo.toml to CRLF in the working copy. autocrlf normalizes it in git, and later edits use a newline-preserving helper (scratch `edit.py`).

## Commits

- 009f64a feat(bench): timing feature, --breakdown and tools/breakdown.py; docs/perf baseline (ticket 02)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
