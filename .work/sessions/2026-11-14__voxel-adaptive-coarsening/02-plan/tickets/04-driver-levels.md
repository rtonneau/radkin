# Ticket 04: driver-levels

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `VoxelRun` holds a current coarsening level, with the state, the `Diffusion` operator and the carried steps on that level's grid.
- [ ] After each accepted splitting step in `advance_to`, if level < max_coarsen_levels and `can_coarsen` passes, the grid drops one level (restrict state, rebuild the diffusion for the coarse grid, reset per-voxel carried steps for the new count); the dt proposal is kept.
- [ ] `inject` refines to level 0 (prolong) before adding the deposit; the Rannacher start behaves as before.
- [ ] With `max_coarsen_levels = 0` every existing driver test passes unchanged.
- [ ] New tests: a homogenising run coarsens and keeps the domain totals within rtol 1e-6 of the same run with coarsening off; an injection after coarsening refines to level 0 and conserves the deposit.

**Files to Touch:**
- `src/voxel/driver.rs`
- `src/voxel/sim.rs`
- `src/voxel/reaction.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::
```

Expected:
All voxel tests pass, the new coarsening tests included.

**Notes:**

The background vector is a single uniform state and does not depend on the grid, so it carries over unchanged. Domain totals must use the current level's voxel volume. Expose the current level (for example `VoxelRun::level()`) for tests and for ticket 5. Count level changes in the run stats, so ticket 6 can report them.
