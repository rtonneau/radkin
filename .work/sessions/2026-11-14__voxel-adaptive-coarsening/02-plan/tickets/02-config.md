# Ticket 02: config

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] `VoxelConfig` parses `max_coarsen_levels` (integer >= 0, default 0) and `coarsen_tol_factor` (float > 0, default 0.1) from `[voxel]`.
- [ ] Config load fails with a `RadkinError::Config` naming the axis when a dim is not divisible by 2^max_coarsen_levels.
- [ ] Unit tests cover defaults, valid values, invalid values and the divisibility error.

**Files to Touch:**
- `src/voxel/config.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::config
```

Expected:
All config tests pass, the new ones included.

**Notes:**

Follow the existing parsing style in `VoxelConfig::from_table`. Check where the grid dims are parsed, so the validation sits next to them. Include the new keys in any config summary or hash that checkpoints compare on resume, if one covers `[voxel]` settings.
