# Ticket 03: fuse-half-reactions

**Status:** ✅ Done

## Local Test Result

Ran `cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo test --release --no-default-features`: all exit 0 (57 lib tests plus integration tests pass). Output in `.scratch/tests/2026-11-14__reduce-computational-burden/evidence/03-checks.txt`.
Bench A/B (old = ticket-02 build, alternating; `evidence/03-ab-16.txt`, `03-ab-1.txt`): 16 threads median 0.512 -> 0.434 s (-15.1%); 1 thread 2.890 -> 2.518 s (-12.9%). Rsoda solves 183205 -> 152510. Gate: max scaled error 0.193 -> 0.193, max total rel error 1.4963e-3 -> 1.4963e-3.

## Review Notes

- AC1: new `strang_pair` in `src/voxel/driver.rs` does R(dt/4) D(dt/2, damped only if the attempt is) R(dt/2) D(dt/2, undamped) R(dt/4); 5 solves per small path, counters confirm (voxel half-steps 211968 -> 176640). Background path and carry use the same `half_reaction` calls, so carried steps stay committed only on accepted steps; the 1- vs 4-thread determinism test and background-vs-full test pass.
- AC2: tests green, gate unchanged.
- AC3: kept, well above 5%. Documented as Lever 2 in `docs/perf/voxel.md`.
Base commit: 5cf38676887e50573c6397d387b3519e1d7e2e98.

## Blockers / Challenges

None

## Commits

- 309ffe7 perf: fuse adjacent voxel half-reactions in step-doubling (ticket 03)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
