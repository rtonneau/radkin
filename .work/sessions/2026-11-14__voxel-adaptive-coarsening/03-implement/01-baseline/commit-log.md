# Ticket 01: baseline

**Status:** ✅ Done

## Local Test Result

`cargo run --release --example voxel_bench` (default case): completes, 0.463 s median, 46 steps, 176640 half-steps, 152510 solves, accuracy lines unchanged (max scaled error 0.193, HO3 total 1.5e-3).
Baseline case: `cargo run --release --example voxel_bench -- --repeat 5 --t-end 1e-5`: median 0.974 s (min 0.908), 57 steps / 0 rejected, 218880 half-steps (24360 background), 194805 solves, nfev 8904705, nst 4437432, max spread 7.2e-6; 1 thread (`RAYON_NUM_THREADS=1`, repeat 3): 5.635 s. Totals match those in docs/perf/voxel.md (deterministic work counters repeated identically across runs).
`cargo fmt --check` and `cargo clippy --release --all-targets -- -D warnings` clean. Raw output in `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/evidence/`. Full `cargo test` not run (no library change).

## Review Notes

Base commit 473dfa6f826e48cb0b6fe8c81cd7e8b7545b21d6, clean tree. Review rounds: 0.
AC1: data/voxel_example.toml (8x8x12, divisible by 4) with `--t-end 1e-5` homogenises (max spread 7.2e-6; at the default 1e-7 s it does not, spread 1.1). No new case file needed. voxel_bench now prints final species totals, per-species spread and a `VOXEL_TOTALS` line.
AC2: wall time, steps, solves, nfev/nst, totals and the exact commands are in the new "Coarsening baseline" section of docs/perf/voxel.md. Wall times are machine-noisy (min 0.908 to 1.016 s over 5 repeats); counters and totals are deterministic.

## Blockers / Challenges

None. Note: the default 1e-7 s case does not homogenise, so the baseline uses `--t-end 1e-5`; the stored accuracy reference does not apply to that t_end.

## Commits

- 0003f80 perf: Coarsening baseline for the voxel bench (ticket 01)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
