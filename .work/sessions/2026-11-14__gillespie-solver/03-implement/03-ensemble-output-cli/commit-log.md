# Ticket 03: ensemble-output-cli

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: 18 unit tests and 13 integration tests (`ssa_exact`, `ssa_beam`, `ssa_cli`) pass, 0 failed.

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 165 passed, 0 failed. `cargo test --release --no-default-features ssa` (serial build, no rayon): 32 passed, 0 failed.

## Review Notes

Base: b54ef3c2fb7ae10bcf5a534969e684a4223b09de. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- `src/ssa/config.rs`: top-level `solver` and the `[ssa]` table (`volume`, `trajectories`, `seed`, `tau_leap`); unknown keys and invalid values are config errors naming the key.
- `radkin run --solver lsoda|ssa` overrides the key (`solver_of`); without key and flag the LSODA code path is the one that ran before (the cell setup was only moved into `build_cell`). The existing CLI tests still pass.
- `src/ssa/ensemble.rs`: trajectories run through `map_ordered` (rayon with the `parallel` feature, serial without) in batches of 256, and the statistics are exact integer sums accumulated in trajectory order. `ssa_cli_output_does_not_depend_on_the_thread_count` compares `concentrations.csv` byte for byte between `RAYON_NUM_THREADS=1` and `4`.
- Output times: `--t-eval`, else 200 log-spaced times (`thinned_times`).
- `concentrations.csv`: LSODA columns first (checked against the header of an LSODA run of the same config), then `<species>_std`; `run.json` records solver, volume, trajectories, seed, tau_leap, n_events.
- LSODA-only options are rejected with their names in the message.
- README section `Stochastic solver (SSA)` added.

`tau_leap = true` is parsed but refused with a message until ticket 5 implements it. A config without beam is allowed with the SSA (reactions only), unlike LSODA.

## Blockers / Challenges

None.

## Commits

- 7282381 feat: SSA ensemble, mean/std outputs and solver selection on radkin run (ticket 03)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
