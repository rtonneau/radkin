# Ticket 02: lean-voxel-solve

**Status:** ✅ Done

## Local Test Result

Ran the Verification Step pieces in order (evidence in
`.scratch/tests/2026-11-14__reduce-computational-burden/evidence/02-*`):

```text
cargo fmt --check                                     -> ok
cargo clippy --release --all-targets -- -D warnings   -> ok (also with --no-default-features)
cargo test --release                                  -> 17 suites ok, 109 passed, 0 failed (02-cargo-test.txt)
cargo test --release --no-default-features --lib -- voxel solver -> 41 passed (serial path)
cargo run --release --example voxel_bench -- --repeat 5   (02-bench-parallel.txt)
VOXEL_BENCH median_s=0.501708 min_s=0.497455 repeats=5 threads=16 voxels=768 steps=46 rejected=0
  voxel_steps=211968 background_voxel_steps=29039 solves=183205 nfev=4219543 nst=2174060
  max_total_err=1.496306e-3 max_scaled_err=1.928750e-1
```

Keep-rule A/B, alternating the ticket-01 build (99e8d37, `02-voxel_bench_base.exe`) with
this build (`02-voxel_bench_final.exe`), medians of the per-run medians
(02-ab-final-16t.txt: 6 x `--repeat 5`; 02-ab-final-1t.txt: 3 x `--repeat 3`):

| threads | before [s] | after [s] | change |
|---|---|---|---|
| 16 | 0.541 | 0.505 | -6.7% |
| 1 (`RAYON_NUM_THREADS=1`) | 3.410 | 2.885 | -15.4% |

Accuracy gate: max per-voxel scaled error 0.1929 (baseline 0.1929), max domain-total
relative error 1.496e-3 HO3 (baseline 1.496e-3). Both >= 5% faster and no worse: kept.

## Review Notes

Base commit: 99e8d37d798a4fdc59256eac8106085d6c80e6a6 (no pre-existing changes).

Acceptance Criteria:
- Dedicated integrate-to-t path: `solver::integrate_to` (no `t_eval`/`Solution`,
  writes `y` in place, returns counters and `h_last`); `SolveScratch` holds the kernel
  `Workspace` and sparse pattern, one per rayon pool thread per half-step (slot by
  `current_thread_index`); per-voxel carried step in `CarriedSteps` on `VoxelRun`
  (follows the accepted trajectory, reset on `inject`, not in checkpoints), used as
  `0.01 * h_last` for the next `first_step`. rsoda's own State/LU buffers are still
  allocated per solve (see Blockers).
- Bit-identity: new tests `integrate_to_matches_solve_with_a_reused_scratch`
  (dense+sparse, reused scratch, bit-equal to `solve` incl. nfev/njev/nst),
  `carried_steps_follow_each_voxel` (parallel step vs per-voxel serial reference,
  background voxels included), `carried_steps_do_not_depend_on_thread_count` and
  `seeded_run_does_not_depend_on_thread_count` (1- vs 4-thread pools, bit-equal state,
  background, steps, carried steps). Existing voxel tests unchanged and green.
- Bench: kept (-6.7% 16 threads, -15.4% 1 thread, gate unchanged), recorded in
  `docs/perf/voxel.md` ("Lever 1").

Checked also: unit-cell `solve` behaviour unchanged (only internal refactor:
`ChemSystem` borrows the workspace; golden tests green); checkpoint format untouched
(`tests/voxel_checkpoint.rs` green); `Workspace` reuse is safe because stale switch
values `sw` are only read for reactions with negative products, which then always
refresh them (`kernel::derive`). njev from ticket 01: not a counter bug, the Jacobian
is never evaluated (no solve leaves Adams); now counted in `ReactionStats::njev`.
Public API changes inside the crate: `reaction_step_background` takes
`&mut CarriedSteps`, `integrate_one` takes `&mut f64` carried step (only the driver
uses them); `reaction_step` keeps its signature (no carry).

## Blockers / Challenges

- First attempt with rayon `try_for_each_init` gave no gain at 16 threads: a counter
  showed 154224 scratch builds for 183205 solves (one per split job). Replaced by one
  slot per pool thread.
- The carried step taken as is (`h_last`) raised nfev (4.65M vs 4.39M): rsoda restarts
  at order 1 and fails the error test with it. Scanned the factor (0.003 to 1, and
  rsoda's own `h0 = 0`); 0.01 gives the fewest RHS calls. The carry lever alone is
  small (-4% nfev); most of the gain is the lean path.
- Reusing rsoda's own buffers (State arrays, dense n x n iteration matrix that is
  never used here) needs an rsoda API. Not needed for the keep rule, so not BLOCKED;
  measured as a follow-up candidate with a local, uncommitted prototype (rsoda copied
  to the session scratchpad, `Lsoda::restart(t0, y0, h0)` resetting `State` in place
  and keeping buffers/linear solver; temporary `[patch]` removed afterwards, rsoda
  checkout untouched): bit-identical in the unit tests, 0.495 s at 16 threads (-12.3%
  vs ticket 01, about -7% on top of this lever) and 2.94 s at 1 thread (-20.6%).
  Numbers in `docs/perf/voxel.md`.
- Timing noise at 16 threads is about +-3% between runs, so the A/B alternates builds
  and uses several runs.

## Commits

- 3fed416 perf: lean per-voxel reaction solve with per-thread scratch and carried step (ticket 02)

## Time Spent

35m (ticket-start.js to ticket-complete.js)
