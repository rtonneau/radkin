# Ticket 04: embedded-error-estimate

**Status:** ✅ Done

## Local Test Result

Final tree (exploration code reverted, only `docs/perf/voxel.md` changed),
evidence in `.scratch/tests/2026-11-14__reduce-computational-burden/evidence/`:

- `cargo fmt --check`: OK. `cargo clippy --release --all-targets -- -D warnings`:
  clean. `cargo test --release`: all green (lib 57, including
  `error_decreases_with_tolerance` and `seeded_run_does_not_depend_on_thread_count`;
  `voxel_checkpoint` 3). `cargo test --release --no-default-features`: all green
  (lib 55). See `04-checks.txt`.
- `cargo run --release --example voxel_bench -- --repeat 5` (`04-final-bench.txt`):
  `VOXEL_BENCH median_s=0.438588 min_s=0.415582 repeats=5 threads=16 voxels=768 steps=46 rejected=0 voxel_steps=176640 background_voxel_steps=24360 solves=152510 nfev=3693056 nst=1896047 max_total_err=1.496305e-3 max_scaled_err=1.928750e-1`.
  This is identical to ticket 03 (the lever was dropped).

Candidate A/B (ticket-03 build `4a46fd0` = `04-voxel_bench_old.exe`, alternating
with the exploration build `04-voxel_bench_explore.exe`; medians of per-run medians;
`04-ab-16.txt`, 4 runs of `--repeat 5`; `04-ab-1.txt`, 3 runs of `--repeat 3`):

| threads | step-doubling | reversed pair, average, K=1 | reversed pair, D/2 R D/2, K=2/3 |
|---|---|---|---|
| 16 | 0.428 s | 0.430 s (+0.4%) | 0.381 s (-10.9%) |
| 1 | 2.515 s | 2.317 s (-7.9%) | 2.101 s (-16.5%) |

The full voxel test suite also passed with each of the two candidates switched on
(`04-candidate-tests.txt`).

## Review Notes

Base commit: `4a46fd04e97d312b10e02c5b1dbac51a98988c42`.

Outcome: **dropped**. `src/voxel/driver.rs` is reverted to the base (step-doubling
and `strang_pair` unchanged). There is no ADR (only for a kept lever) and no README
change. The numbers and the analysis are in `docs/perf/voxel.md`, section "Lever 3".

Against the Acceptance Criteria:

- AC1 (replace step-doubling): implemented and measured in an exploration build,
  not kept. Candidates: Strang vs Lie sharing the first half-reaction (exponent -1/2)
  and Strang vs reversed Strang `R/2 D R/2` / `D/2 R D/2` (exponent -1/3), keeping
  S1, S2 or their average, with an error constant K scanned (`04-scan-errk.txt`).
  The reversed pair costs about 1.6 Strang-equivalents per attempt (one solve over
  dt = 1.22 solves over dt/2).
- AC2 (ADR): not created, because the lever was not kept.
- AC3 (tests, gate): green on the final tree, and the gate is unchanged. Both
  candidates that pass the official gate also pass the full test suite.
- AC4 (bench >= 5% or dropped and recorded): dropped and recorded.
  - Strang vs Lie: 135 steps, 0.722 s vs 0.455 s.
  - Average, K=1: passes the gate, but +0.4% at 16 threads.
  - D/2 R D/2, K=2/3: formally passes the gate (1.210e-3 / 0.154) and the keep
    rule (-10.9% / -16.5%). But K was tuned on the gate case.

Why that last candidate is dropped anyway: held-out end times of the same case
(new references `04-voxel_reference_t1e-8.json` and `_t1e-6.json`, written with
step-doubling at rtol/atol 1e-6; `04-scan-tend*.txt`).

- At 1e-8 it gives 2.685e-3 / 0.533, against step-doubling's 1.237e-3 / 0.232.
- At 1e-6 it gives 1.021e-3 / 0.050, against 4.093e-4 / 0.027.
- Matching step-doubling needs K of about 3, which costs more nfev than
  step-doubling.

Root cause (`04-sd-keep.txt`): at step-doubling's own steps, keeping one Strang
step of dt gives about 4x the error of the kept two-half-step solution. So any
embedded pair needs about 2x the steps. The reversed-pair difference
`(P - Q)/8 dt^3` is not proportional to either ordering's error and can cancel.
At 1e-7 only, a single D/2 R D/2 step is unusually accurate.

Not tried, recorded as a possible follow-up: step-doubling with the first
half-reaction of both paths taken from one solve (tcrit stop at t + dt/4), 4
solves instead of 5.

## Blockers / Challenges

- No candidate gave a robust gain. The one that formally passed the gate and the
  keep rule (D/2 R D/2, K=2/3) was 1.8-2.5x less accurate than step-doubling at
  t_end 1e-8 and 1e-6. I dropped it as overfit to the gate case rather than
  keeping it on the letter of the rule. Reproduce any candidate with the
  exploration build: `RADKIN_SPLIT=rev_s2|rev_avg|rev|lie RADKIN_ERRK=<K>
  04-voxel_bench_explore.exe`.
- `evidence/04-explore-timing-16t.txt` is a preliminary same-binary timing. Its
  `steps=` column is mislabeled (it shows background voxel steps). It is annotated
  in the file and superseded by `04-ab-16.txt` / `04-ab-1.txt`.

## Commits

- fc39c67 docs: embedded voxel error estimate measured and dropped (ticket 04)

## Time Spent

20m (ticket-start.js to ticket-complete.js)
