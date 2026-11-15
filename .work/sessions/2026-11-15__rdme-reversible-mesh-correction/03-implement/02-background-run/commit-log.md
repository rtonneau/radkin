# Ticket 02: background-run

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`: exit 0, clippy clean, all 25 test binaries ok (lib 146 passed plus the new `parses_background_run` and `ssa_rdme_background_run_nets_a_seedless_run_to_zero`). Log: `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/evidence/02-verify.log`.

## Review Notes

- AC1: `[ssa] background_run` is an `Option<bool>` (default on); non-bool values are rejected; any explicit value on a unit-cell run is an error in `require_volume` (test `parses_background_run`).
- AC2: `run_rdme` runs a second ensemble on the same model and record times, with no seeds or snapshots, seed `cfg.seed ^ 0x9E3779B97F4A7C15` (`background_seed`). (Superseded by the correction pass below: only G is background-net.) Snapshot fields are unchanged. No seed field is kept for the background.
- AC3: `snapshots.json` gets `background_run` (seed, per-species noise = SEM of the background run's mean change over the whole run, event counts); it is null when opted out.
- AC4: `ssa_rdme_background_run_nets_a_seedless_run_to_zero` (seedless 16^3 x 20 nm water box, 200 trajectories): all net means within 4 SEM of zero; opt-out and JSON output checked.
- First pass only (reverted): the existing gross-total comparisons against LSODA (`ssa()` helper in `tests/ssa_rdme.rs`, and `tests/ssa_rdme_aerated.rs`, outside the listed files) now set `background_run = Some(false)`, since the default would subtract the background from totals.

Correction pass 1 of 1 (parent): N columns and gross totals are kept as before; the background run's totals are subtracted only for the net (G) values (`net_mean` = totals minus background totals, `net_std` combined in quadrature), as in `src/voxel/output.rs`. The `background_run = Some(false)` opt-outs in `tests/ssa_rdme.rs` (`ssa()` helper) and `tests/ssa_rdme_aerated.rs` were removed; no test needs one (the only remaining `Some(false)` is the opt-out check in the new test). The seedless net-zero test now reads `net_mean`/`net_std` only. A unit-cell run now errors only on `background_run = true`; `false` is accepted (config test updated). Rerun of the Verification Step: see Local Test Result.

## Blockers / Challenges

None remaining (the first pass's N-column subtraction was reverted in the correction pass).

### Parent review and acceptance

- Base 6454823b8e6537294f25f005ffb2e5a83e68e997; pre-existing changes: gps session state only. Correction passes used: 1 (parent: N_ columns and totals stay gross and only the G numerator is net, as in src/voxel/output.rs; a unit cell errors only on `background_run = true`).
- Independent review (fresh reviewer, base 6454823): no actionable findings. P3 `docs/ssa.md` not updated: that is ticket 03.
- Accepted: every Acceptance Criterion is met, and the full Verification Step passes (`evidence/02-verify.log`). Requested model sonnet/medium.

## Commits

- cb0fb5c feat: background-net RDME yields from a seedless background run (ticket 02)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
