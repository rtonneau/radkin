# Ticket 03: coarsen-primitives

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] A new `src/voxel/coarsen.rs` provides: the coarse `VoxelGrid` of a grid (dims/2, spacing*2, same origin), `restrict` (volume average of each 2x2x2 block, every species), `prolong` (piecewise constant) and `can_coarsen` (per-species spread over the 8 children and their face neighbours below `factor * (atol + rtol*|mean|)` for every block).
- [ ] Unit tests: restrict conserves the domain totals of every species to round-off (amount = concentration * voxel volume); prolong(restrict(x)) is exact for blockwise-constant x; restrict(prolong(y)) == y; a uniform field can coarsen; a single spike cannot.
- [ ] Primitives do not allocate per voxel (buffers passed in or allocated once per call).

**Files to Touch:**
- `src/voxel/coarsen.rs`
- `src/voxel/mod.rs`
- `src/voxel/grid.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::coarsen
```

Expected:
All coarsen tests pass.

**Notes:**

`VoxelState` stores species-major fields (see `field`, `gather`, `scatter` in `grid.rs`). Concentrations are per volume, so the volume average of equal-size children is the plain mean. Face neighbours at the domain boundary: skip them (zero-flux walls). Keep this module independent of the driver.
