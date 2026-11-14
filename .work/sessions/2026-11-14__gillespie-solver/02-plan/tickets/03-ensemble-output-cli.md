# Ticket 03: ensemble-output-cli

**Acceptance Criteria:**
- [ ] `[ssa]` table: `volume` (litres, required for the unit cell), `trajectories`, `seed`, `tau_leap` (default false); top-level `solver = "lsoda" | "ssa"` (default lsoda); invalid values are config errors with a clear message
- [ ] `radkin run --solver <lsoda|ssa>` overrides the config key; without `solver` and without the flag the LSODA run is unchanged
- [ ] The ensemble runs N trajectories in parallel (rayon, `parallel` feature) and reduces them in trajectory order: the output is byte-identical with 1 thread and with many
- [ ] Output times: `--t-eval` when given, else log-spaced points over the span
- [ ] `concentrations.csv` has the LSODA layout (t, dynamic then fixed species, µmol/L, ensemble mean) followed by one `<species>_std` column per dynamic species; `run.json` records solver, volume, trajectories, seed, tau_leap and the event count
- [ ] LSODA-only flags that make no sense with SSA (`--mode chunked`, `--production`, ...) are rejected with a message rather than ignored
- [ ] README documents the solver selection and the `[ssa]` table
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/mod.rs`
- `src/ssa/config.rs`
- `src/ssa/ensemble.rs`
- `src/ssa/output.rs`
- `src/bin/radkin.rs`
- `tests/ssa_cli.rs`
- `README.md`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including the CLI test that runs `radkin run --solver ssa` twice with different thread counts and compares the files.

**Notes:**

`Config::from_table` ignores unknown root keys, so `solver` and `[ssa]` are read from the root TOML file like `[voxel]`. Thread count for the test: `RAYON_NUM_THREADS`.
