# Ticket 02: beam-sources

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: 13 unit tests and 9 integration tests (`ssa_exact`, `ssa_beam`) pass, 0 failed.

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release` all passed (count in the completion message).

## Review Notes

Base: 372b9fdc4968830eda4b290707314c5993ff74f9. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- Sources: `SsaModel::sources` holds `yield * N_A V` per gray for every radiolytic reaction of a dynamic species (built in ticket 1's model); the stepper multiplies by the segment dose rate. `ssa_sources_are_built_from_the_g_values` checks the value against 2.6 / (100 e) per joule.
- Piecewise-constant dose rate: `src/ssa/schedule.rs::beam_segments` cuts the span at the rows of `Beam::make_matrix` and takes each row's dose rate; `run_trajectory` stops at a segment end, discards the pending event and redraws (memoryless, exact).
- Continuous and pulsed beams: ensemble mean within 4 SEM of `yield * dose * N_A V` at several times, variance equal to the mean (Poisson) within its sampling error, nothing produced between pulses.
- No beam: zero events.

`model.rs` and `exact.rs` needed no change in this ticket: the source handling was written with the stepper core in ticket 1.

## Blockers / Challenges

None.

## Commits

- f61b2a0 feat: SSA beam sources as zero-order reactions on a piecewise-constant dose rate (ticket 02)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
