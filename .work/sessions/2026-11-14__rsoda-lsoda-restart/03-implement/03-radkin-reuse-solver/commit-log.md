# Ticket 03: radkin-reuse-solver

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`:
exit 0, fmt clean, clippy clean, 134 tests passed / 0 failed (17 test binaries, golden tests included,
tolerances unchanged). Full output: `.scratch/tests/2026-11-14__rsoda-lsoda-restart/evidence/03-verification.txt`.

Extra checks: `cargo clippy --release --no-default-features --all-targets -- -D warnings` and
`--features timing` clean; `cargo test --release --no-default-features` 132 passed / 0 failed.

Voxel results: all 67 output files of `radkin voxel data/voxel_example.toml` (checkpoint, snapshots,
yields) are byte-identical between the pre-change build (a91c033, rsoda f5d1312), the rev bump alone and
the final build, at 16 threads and at 1 thread (`evidence/03-voxel-bitwise.txt`). All 132 voxel_bench runs
give the same work and accuracy (46/0 steps, 152510 solves, nfev 3693056, nst 1896047, max total error
1.496305e-3, max scaled error 1.928750e-1).

Bench (`evidence/03-ab-summary.txt`, raw `03-ab-*-raw*.txt`): after vs rev bump alone -2.1/-2.4/-4.7% at
16 threads, -6.0/-8.0/-9.7% at 1 thread; against a same-code control using `Lsoda::new` -4.8% (16) and
-7.9% (1). Against the pre-session build: +0.1% to +2.2%. Keep rule applied: kept (recorded in
`docs/perf/voxel.md`, "Solver restart").

## Review Notes

- AC1: `Cargo.toml` pins rsoda `rev = "2255091"` (full sha 22550912c25a..., on origin/feat/lsoda-restart);
  `Cargo.lock` updated; no `[patch]`.
- AC2: `solver::Integrator` owns a `ChemSystem` (kernel `SolveScratch`) and one `Lsoda` built by the same
  `new_lsoda` (same tolerance, Jacobian kind, hmax); `integrate_to` calls `Lsoda::restart(t0, y, first_step)`,
  so radkin's `first_step` (CARRY_FACTOR x carried step, or `opts.first_step`, `h / 1e3` cap) goes in as `h0`.
  In `reaction_step_background` each pool thread builds one on its first voxel (slot per
  `current_thread_index`, as the scratch before) and restarts it per voxel; the background solve takes a
  solver already built (or builds one if no voxel was active). `restart_warm` is not used. `solve`
  (run/chunked/single-pulse drivers) and `integrate_to` keep `Lsoda::new`; `integrate_one` (one solve) too,
  so the carried-steps test compares the restart path against `new()` bit for bit.
- `ChemSystem` is generic over `K: BorrowMut<SolveScratch>` (borrowed for `integrate_to`, owned for
  `Integrator`), so no `Arc`, no self-reference; the kernel still does not allocate. The stepping loop
  moved to a shared `step_to`.
- New unit-test coverage: `integrate_to_matches_solve_with_a_reused_scratch` also restarts one
  `Integrator` over states x1, x3, x1 in dense and sparse mode and asserts equal states and counters.
- AC3/AC4: see Local Test Result. Keep decision: the rule asks >= 5% on the targeted cases with none
  slower; one thread passes in every session (and against the same-code control), 16 threads gains 2-5%,
  never slower. Reviewer note: against the pre-session build there is no net gain, because the rev-bump-only
  build runs 2-13% slower than it (layout difference, see Blockers). If the reviewer reads the rule as
  "16 threads must gain 5%", the revert is `git checkout -- src/solver.rs src/voxel/reaction.rs` plus
  rewording the voxel.md section; the rev bump can stay.

### Parent review record

- Base: a91c033300abb3ca0b4281536ef161f0043cd5f9. Pre-existing uncommitted files: only gps session state under .work/sessions/.
- Independent reviewer (fresh context, opus/high): REVIEWED, 3 findings, no correctness defects (options identical to new_lsoda, failed solve does not poison later voxels, per-thread slots sound, other drivers unchanged, no [patch], kernel allocation-free).
  - P2 docs/perf/voxel.md keep claim within noise (16 threads below 5%, 1-thread margin from one control session). Disposition: section reworded to 'Kept provisionally'; the final keep or revert is decided in ticket 04 with repeated A/B sessions (new vs restart vs restart+warm). If restart alone fails there, the code is reverted in a follow-up commit and the rev bump stays.
  - P3 voxel.md 0.06 us figure: fixed (0.07-0.1 us).
  - P3 solver.rs Integrator::new doc on opts.first_step: fixed.
- Gaps accepted: no radkin unit test of reuse after a failed solve (covered by rsoda tests/restart.rs and per-call drop); background reuse checked end-to-end only.
- Parent re-ran fmt, clippy -D warnings, cargo test --release: 134 passed, 0 failed.
- Accepted with the provisional keep above.

## Blockers / Challenges

- The rev bump alone made the voxel bench 2-6% (16 threads) / 7-13% (1 thread) slower than the pre-session
  build, with identical results and counts. Checked it is not an rsoda regression: rsoda 2255091 only changes
  setup code (microbench `Lsoda::new`+drop 0.68 -> 0.75 us; `restart` 0.11 us), and a run-mode control
  (examples/bench.rs, Labarbe 1e-7..10 s) is 5.6% faster on the new rev. So it is a build-to-build code layout
  difference; the restart gain cancels it, which is why the final build equals the pre-session one. To
  separate the restart from layout, a control build (`03-voxel_bench_afternew.exe`: final code with
  `Lsoda::new` instead of `restart`, temporary edit reverted afterwards, not committed) was measured: it
  matches the rev-bump build, and restart is -4.8% / -7.9% against it.
- The expected 16-thread allocator-contention gain (prototype -7%) did not appear beyond the per-solve
  saving; the 1-thread gain is larger than the prototype's (none).
- A first 16-thread A/B run lost its output file (wrong relative path); it was rerun.

## Commits

- 346f8ad perf: Reuse one rsoda solver per thread in the voxel reaction step (ticket 03)

## Time Spent

34m (ticket-start.js to ticket-complete.js)
