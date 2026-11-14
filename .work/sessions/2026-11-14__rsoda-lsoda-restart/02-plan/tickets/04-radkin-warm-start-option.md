# Ticket 04: radkin-warm-start-option

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] A `[voxel]` config option (off by default) turns on rsoda's warm start in the voxel reaction step; the first step still comes from `CARRY_FACTOR` times the carried step.
- [ ] A test covers parsing the option and its default.
- [ ] The voxel bench with the option on and off reports nfev and wall time, and the accuracy against the reference run, in `docs/perf/voxel.md`. If it does not need fewer RHS calls and less wall time at equal accuracy, the option is removed and the result is recorded.
- [ ] radkin checks are clean.

**Files to Touch:**
- `src/voxel/config.rs`
- `src/voxel/reaction.rs`
- `README.md`
- `docs/perf/voxel.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All clean; all tests pass.

**Notes:**

If ticket 2 pushed new rsoda commits after ticket 3's pin, bump the `rev` again here.
