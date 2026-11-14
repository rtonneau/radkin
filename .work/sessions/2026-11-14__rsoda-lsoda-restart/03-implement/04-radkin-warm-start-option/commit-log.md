# Ticket 04: radkin-warm-start-option

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` on the
final tree: exit 0, fmt clean, clippy clean, 134 passed / 0 failed (17 test binaries, golden tests at unchanged
tolerances). Output: `.scratch/tests/2026-11-14__rsoda-lsoda-restart/evidence/04-verification.txt`.
`src/solver.rs` and `src/voxel/reaction.rs` equal `a91c033` (`git diff --quiet a91c033`); Cargo.toml/lock keep the
rsoda `2255091` pin. Final voxel_bench work identical to before: 46/0 steps, 152510 solves, nfev 3693056, nst
1896047, max total err 1.496305e-3, max scaled err 1.928750e-1.

The warm-start option was implemented and tested before the measurement (`[voxel] warm_start`, bool, default
false; `Integrator::with_warm_start`; plumbed through `VoxelRun::warm_start` and `sim.rs`; `voxel_bench
--warm-start`; config test for parsing, default and type errors) and passed build; it was then removed by the
verdict below.

## Review Notes

- AC1 (option off by default, first step CARRY_FACTOR x carried step): implemented, measured, then removed
  per the ticket's own rule (AC3).
- AC2 (parse + default test): written (`warm_start_is_off_by_default_and_parsed`), removed with the option.
- AC3 (bench on/off, nfev, wall, accuracy in docs/perf/voxel.md): done. Two bench methods, 3 sessions each at
  1 and 16 threads, orders new/restart/warm, warm/restart/new, restart/new/warm, short case (t_end 1e-7) and late
  case (t_end 1e-5, own reference written to `.scratch/.../04-voxel_reference_1e-5.json`).
  - Warm vs restart (same binary): short case, no solve ends on BDF, so the work is identical and the wall-time
    difference is noise (-0.8% to +3.2%). Late case: nfev -6.3% to -7.0% (16t) / -10.6% (1t) at equal accuracy
    (2.556072e-4 vs 2.556073e-4), wall -1.2% to -1.5% (16t), -3.3% to +0.2% (1t). Also nondeterministic at 16 threads
    (12 distinct nfev values in 12 runs: the method carried is that of the thread's previous voxel). Verdict: REMOVED.
  - Restart vs new: separate builds (control = temporary edit) gave restart 0.5-11% SLOWER, the opposite sign of
    ticket 03's control (5-8% faster), so cross-build layout noise exceeds the effect. Decisive run in one binary
    (env var switch, temporary, not committed): restart -4.5/-4.6/-5.0% (16t short, pooled -4.4%), -2.6..-3.4%
    (16t late), -6.1..-9.6% (1t short), -1.5..-4.8% (1t late); never slower. Keep rule (>= 5% on targeted cases)
    fails at 16 threads on the default case and on the late case. Verdict: REVERTED (`src/solver.rs`,
    `src/voxel/reaction.rs` back to a91c033; rev bump kept).
  - docs/perf/voxel.md: the section is renamed "- reverted" (the Lever 1 link anchor was updated), the
    "Kept provisionally" wording is replaced with both verdicts, tables, and the reason the two-binary method was
    not used for the decision.
- AC4: checks clean (see Local Test Result).
- README.md: not changed (the option it would have documented was removed).
- Reviewer note: the restart verdict is a strict reading of the 5% rule. In one binary restart was never slower
  and saves about 1.6 us per solve; if the rule is read as "1 thread is the target", restoring ticket 03 is
  `git checkout 346f8ad -- src/solver.rs src/voxel/reaction.rs` plus rewording the voxel.md section.

### Parent review record

- Base: 54250b8e179561d8a25fec5ce41711e0fa3b0ab9; pre-existing uncommitted files: gps session state only.
- Independent review skipped: the code diff is a revert of src/solver.rs and src/voxel/reaction.rs to a91c033 (git diff a91c033 -- src is empty; only the Cargo.toml rev pin differs), plus docs/perf/voxel.md. The parent checked that diff and re-ran fmt, clippy -D warnings, cargo test --release: 134 passed, 0 failed.
- Verdicts accepted per the keep rule: warm start removed (no wall-time gain, nondeterministic at 16 threads); restart reverted (16-thread default case -4.4% to -5.0%, under 5%). The rsoda pin at 2255091 stays (restart and restart_warm remain available in rsoda).
- Note: the warm-start option was built and tested, then removed under the ticket's own rule, so README.md and config.rs are unchanged.

## Blockers / Challenges

- Build-to-build layout noise (about +/-10%) is larger than the restart effect: the two-binary A/B reversed sign
  against ticket 03's. Resolved by a one-binary control (`evidence/04-control-switch.diff`, reverted).
- `restart_warm` cannot be made deterministic per voxel with the current rsoda API (no restart on a given method);
  noted in voxel.md as what a deterministic version would need.
- Temporary edits (control builds, BDF counter) were all reverted; the diffs are kept in evidence/.

## Commits

- ede6e41 perf: Revert the voxel solver restart, drop the warm start (ticket 04)

## Time Spent

32m (ticket-start.js to ticket-complete.js)
