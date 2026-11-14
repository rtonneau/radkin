# Ticket 01: voxel-bench-and-gate

**Status:** ✅ Done

## Local Test Result

`cargo run --release --example voxel_bench -- --repeat 3` (exit 0, after `--write-reference`, which took 5.2 s, 578 steps):

```
time: median 0.5736 s, min 0.5704 s over 3 repeats
steps: 46 accepted, 0 rejected; voxel half-steps 211968 (29039 background); rsoda solves 183205; nfev 4389691, nst 2202873
max domain-total rel. error: 1.496e-3 (HO3); max per-voxel scaled error: 1.929e-1
VOXEL_BENCH median_s=0.573639 ... max_total_err=1.496305e-3 max_scaled_err=1.928750e-1
```

`RAYON_NUM_THREADS=1`: median 3.365 s, same counts and accuracy. Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (105 passed, 0 failed). Outputs kept in `.scratch/tests/2026-11-14__reduce-computational-burden/evidence/` (reference.txt, bench_parallel.txt, bench_1thread.txt).

## Review Notes

Base commit 1465054892b2e47e50c85239be7318b4c5649463, no pre-existing changes. Review rounds: 0.

- AC1: `examples/voxel_bench.rs` runs the case `--repeat N` times, prints median/min time, accepted/rejected steps, voxel half-steps, rsoda solves, nfev, nst, and a `VOXEL_BENCH key=value` line. Solves/nfev/nst are not in `VoxelOutput`, so `src/voxel/reaction.rs` got read-only relaxed atomic counters (`ReactionStats`, `stats()`, `reset_stats()`), incremented in `reaction_step_background`; the integration is unchanged (voxel tests and full suite pass). Counters are process-global, so they assume one voxel run at a time.
- AC2: `--write-reference` runs the case at splitting rtol 1e-6 (atol tightened to 1e-6 too, flags `--reference-rtol/--reference-atol`) and writes the final state to the scratch dir; normal runs print per-species domain-total relative error and max per-voxel scaled error (using the run's own rtol/atol).
- AC3: baseline in `docs/perf/voxel.md`, linked from `docs/perf/README.md`.
- Chose the configuration's t_end = 1e-7 s (0.57 s on 16 threads, 3.4 s on 1).
- Dropped `njev` from the output: it read 0 for all 183k voxel solves, cause not investigated.
- Gap: the reference tightens atol as well as rtol (ticket said rtol only); documented.

## Blockers / Challenges

None (a heredoc write of the example failed in the shell once; rewrote it with the Write tool).

## Commits

- ea03246 perf: voxel benchmark and accuracy gate (ticket 01)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
