# Ticket 01: voxel-grid-config

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] A `VoxelGrid` (dims, spacing, origin, SoA state storage, index/position helpers) exists in `src/voxel/grid.rs`.
- [ ] TOML `[grid]` section parses with validation errors (`RadkinError`) for zero dims or non-positive spacing.
- [ ] Optional per-species `diffusion` coefficient (m²/s) parses; species without one are immobile.
- [ ] Existing configs and golden tests are unchanged.

**Files to Touch:**
- `src/voxel/mod.rs`
- `src/voxel/grid.rs`
- `src/config/mod.rs`
- `src/chemistry/species.rs`
- `src/lib.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::grid && cargo test --release
```

Expected:
All tests pass, including new grid/config tests.

**Notes:**

Default D values follow Geant4-DNA (read the G4 source `G4EmDNAChemistry*` to copy them; do not guess).
