# Ticket 01: tonneau-config

**Status:** ✅ Done

## Local Test Result

`cargo test --release --test golden_runs paper`: 2 passed (paper_points_match_golden, tonneau_2025_config_matches_original_paper_points). `radkin run data/tonneau_2025.toml --mode single-pulse --t-span 0 500`: 3390 points written. `cargo fmt --check` clean. The test went in golden_runs.rs instead of config.rs because the paper-point harness lives there.

## Review Notes

Base 27c3d411d8f855f90b8dbdd06b3228ff3c0b57e4, clean tree. Inline, localized change; no independent reviewer. The new test reproduces the 12 golden "original" points at the golden rtol from the file alone: variable G-values, and LET taken from the config beam. All acceptance criteria are met.

## Blockers / Challenges

The radiopyo beam type key is `singlepulse`, not `single_pulse`. Fixed.

## Commits

- 5bf74ee feat: Bundle data/tonneau_2025.toml, the Tonneau 2025 original-runs inputs (ticket 01)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
