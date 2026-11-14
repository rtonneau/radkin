# Ticket 05: io-resume

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] Snapshots and the CSV/JSON outputs are always on the fine grid (the current state prolonged), whatever the level.
- [ ] Checkpoints store the fine-grid state; the checkpoint format is unchanged (existing checkpoint tests pass).
- [ ] On resume the level is recomputed from the state (coarsen repeatedly while `can_coarsen` passes, up to max_coarsen_levels); a test shows resume from a mid-run checkpoint matches an uninterrupted coarsened run within rtol 1e-6.

**Files to Touch:**
- `src/voxel/sim.rs`
- `src/voxel/checkpoint.rs`
- `src/voxel/output.rs`
- `src/voxel/driver.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::
```

Expected:
All voxel tests pass, checkpoint and resume tests included.

**Notes:**

Prefer a single accessor on `VoxelRun` that returns the fine-grid state, used by snapshots and checkpoints alike, so no caller sees a coarse state.
