# Ticket 02: ho3-diffusion

**Status:** ✅ Done

## Local Test Result

`cargo test --release --lib voxel::config`: 3 passed. Full suite `cargo test --release`: all ok, 0 failed. `cargo clippy --release --all-targets -- -D warnings` and `cargo fmt --check` clean.

## Review Notes

Base fc72dd5 (after ticket 01). Inline, localized change; no independent reviewer. DEFAULT_DIFFUSION has 16 entries, HO3 = 2.0e-9 with a comment saying Geant4-DNA has no HO3 (checked: no HO3 in geant4-v11.4.3 source/processes/electromagnetic/dna). parses_grid_and_diffusion asserts the default. The README line is updated.

## Blockers / Challenges

None.

## Commits

- b023167 feat: Default HO3 diffusion coefficient (O3's 2.0e-9, absent from Geant4-DNA) (ticket 02)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
