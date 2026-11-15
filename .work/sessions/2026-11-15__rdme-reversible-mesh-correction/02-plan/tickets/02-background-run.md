# Ticket 02: background-run

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `[ssa] background_run` (bool, default true for voxel runs) is parsed. Unknown values are rejected, and setting it on a unit-cell run is an error, the same way `rate_correction` is.
- [ ] A voxel SSA run with the flag on runs a seedless trajectory on the same grid and settings, with an RNG seed derived from the run's seed, and subtracts its counts so that `yields.csv` is background-net.
- [ ] The output metadata records the background run and reports its noise.
- [ ] Tests cover the config parsing and show the net yields of a seedless run are zero within noise.

**Files to Touch:**
- `src/ssa/config.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/voxel_output.rs`
- `tests/ssa_rdme.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All tests pass; clippy clean.

**Notes:**

This mirrors how the deterministic voxel path integrates the seedless background (session 2026-11-13). Keep memory in mind: the background run needs no seed field.
