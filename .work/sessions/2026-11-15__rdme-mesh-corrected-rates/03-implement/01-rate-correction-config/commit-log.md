# Ticket 01: rate-correction-config

**Status:** ✅ Done

## Local Test Result

`cargo test --release --lib ssa::config`: 4 passed (parses_solver_and_ssa_table, defaults_to_lsoda, unit_cell_refuses_an_explicit_mesh_correction, rejects_invalid_values). `cargo clippy --release --all-targets -- -D warnings` clean; `cargo fmt --check` clean after `cargo fmt`.

## Review Notes

Base: 1422cb231ec835822573dae141744c2919d87b35. Pre-existing untracked files (not owned): `.work/sessions/.pending-seeds.json`, `.work/sessions/scout-reports/`, `docs/literature/`.

- `RateCorrection::{None, Hellander}` and `SsaConfig::rate_correction: Option<RateCorrection>` (unset = `None`, the default with fallbacks), re-exported from `ssa`.
- Parse accepts "hellander"/"none"; other strings and non-strings are config errors naming `ssa.rate_correction`; the unknown-key message lists it.
- The unit-cell check sits in `require_volume`, which `run_ensemble` calls, so the library and CLI paths both refuse an explicit "hellander".
- Module doc shows the key. Struct literals in `tests/ssa_tau.rs` and `examples/ssa_bench.rs` gained the field.

## Blockers / Challenges

A heredoc-fed Python edit lost a level of backslashes in Git Bash and aborted on its first assertion before writing; redone with Edit.

## Commits

- 370aa4e feat: [ssa] rate_correction key, refused on unit-cell runs (ticket 01)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
