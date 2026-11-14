# Ticket 06: acceptance

**Status:** ✅ Done

## Local Test Result

Verification Step (base commit dab5917, working tree with this ticket's changes):

```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example voxel_bench
```

Result: exit 0. fmt clean, clippy clean, 134 tests passed / 0 failed, and the bare bench ran
the default case (t_end 1e-7, coarsening off): median 0.452 s, 46 steps, 176640 half-steps,
152510 solves, nfev 3693056, accuracy vs reference unchanged (max total error 1.496e-3 (HO3),
max scaled error 0.193). Output: `.scratch/tests/2026-11-14__voxel-adaptive-coarsening/evidence/06-verify.txt`.

The bare command does not coarsen (the example config has `max_coarsen_levels` 0 and the
case never homogenises by 1e-7). The coarsened-vs-fine comparison is the new bench mode,
`--max-coarsen-levels N` (it runs the fine grid alongside, alternating with the coarsened
run, and prints both times, the level history and the accuracy against the fine run):

```bash
cargo run --release --example voxel_bench -- --repeat 5 --t-end 1e-5 --max-coarsen-levels 2
```

Results at 16 threads, `--repeat 5`, factor 0.1 (default): t_end 1e-5 fine 0.885 s, coarsened
0.661 s (x1.34), levels `1e-13:1, 3e-13:2, 1e-12:0, 1.585e-6:1, 2.929e-6:2`, half-steps
218880 -> 197040, nfev 8904705 -> 6171868, max domain-total rel. error vs fine 2.7e-11 (O_minus;
criterion rtol 1e-6), max field difference 1.1e-3 merge tolerances (1.1e-4 control
tolerances). One thread: 5.71 s -> 4.11 s (x1.39). t_end 1e-4: x1.95, total error 1.1e-9.
t_end 1e-7 and 1e-6: x1.01, no coarsening before the end. `coarsen_tol_factor` sweep at 1e-5
(0.01 .. 100) in `docs/perf/voxel.md`; total error 2.7e-11 (0.1), 1.8e-9 (0.3), 4.0e-8 (1),
9.7e-7 (3), 1.7e-6 (10, fails), 2.2e-4 (100); speedup x1.34, 1.44, 1.51, 1.66, 1.73, 2.19.
All raw runs are in the scratch evidence folder (`06-run-*.txt`).

## Review Notes

- Base commit dab5917 (HEAD), tree clean apart from session bookkeeping. Changed files:
  `examples/voxel_bench.rs`, `docs/perf/voxel.md`, `README.md`, `src/voxel/config.rs` (doc
  comment only). No library behaviour changed; no review rounds yet (0).
- AC1 (baseline case runs coarsened, totals within rtol 1e-6, fields within merge tolerance):
  met at factor 0.1 and 1 (2.7e-11 / 4.0e-8 totals; fields 1.1e-3 / 1.3e-4 merge tolerances),
  at 16 and 1 threads and at t_end 1e-5, 1e-4.
- AC2 (factor tuned): default kept at 0.1. Factor 1 also passes (25x margin on the totals,
  x1.51 vs x1.34), but the totals error grows about as factor squared (3 is at the limit, 10
  fails), the margin comes from one case, and the extra 13% is small. Recorded in the docs.
- AC3 (docs/perf): section "Global coarsening" with times, speedup, level history, the
  keep-rule decision (kept, off by default) and the local-AMR assessment: global coarsening
  cannot touch the window before the field homogenises (until 1.6e-6 s here, about 70% of the
  fine cost of a 1e-5 run), so local AMR would be the way to speed that window; not measured.
  Honest summary: the speedup is modest on short runs (x1.01 at t_end 1e-7 and 1e-6, the
  default case), x1.34 at 1e-5 and x1.95 at 1e-4.
- AC4: README documents `max_coarsen_levels` and `coarsen_tol_factor` in the `[voxel]` block
  and a paragraph; full CI checks clean (above).
- Bench change: `--max-coarsen-levels`, `--coarsen-tol-factor`, and in coarsened mode an
  alternating fine run with times, level history, per-species total error and field
  difference; `VOXEL_BENCH` line gains `fine_median_s`, `speedup`, `coarse_levels`,
  `coarse_max_total_err`, `coarse_max_scaled_err` only in that mode, so existing script
  output is unchanged. Alternation was added after a first sweep (coarsened run first) showed
  a warm-up bias (the first run was 10% slower at t_end 1e-7); the first-sweep files are in
  `06-first-sweep/` in the scratch folder and are not cited.
- Reviewer observation confirmed: the 5 level changes at 1e-5 include 2 coarsenings of the
  uniform pre-seed state (1e-13, 3e-13) refined at the seed (1e-12); only 2 changes are real.
  Noted in docs/perf; no solves are spent in that window so no code change.

## Blockers / Challenges

None. Notes: (1) the first sweep had a measurement bias (coarsened run always first in the
process), fixed by alternating the runs and redoing the sweep. (2) A Python heredoc edit wrote
CRLF into `docs/perf/voxel.md` and `examples/voxel_bench.rs` (repo is `eol=lf`); converted back
to LF. (3) The bare Verification command does not exercise coarsening (off by default); the
coarsened comparison needs `--max-coarsen-levels`.

## Commits

- 722428f perf: Measure global voxel coarsening and document its keys (ticket 06)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
