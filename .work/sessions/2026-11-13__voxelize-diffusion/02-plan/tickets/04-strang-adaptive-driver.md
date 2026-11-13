# Ticket 04: strang-adaptive-driver

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] A `VoxelRun` driver does Strang steps R(dt/2) D(dt) R(dt/2).
- [ ] Step-doubling error control with configurable rtol/atol, dt growth/shrink limits, and exact landing on event times (injections, snapshots).
- [ ] Test: a uniform seed through the full driver matches Phase 1 within rtol 1e-6.
- [ ] Test: on a small reaction-diffusion case, the error decreases as tolerance tightens.

**Files to Touch:**
- `src/voxel/driver.rs`
- `src/voxel/mod.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::driver
```

Expected:
Driver tests pass.

**Notes:**

The accepted state is the two-half-step solution.
