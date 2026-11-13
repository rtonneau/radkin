# Ticket 03: voxel-reaction-step

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] The reaction half-step integrates every voxel over [t, t+h] with rsoda using the existing kernel (no beam source).
- [ ] A rayon-parallel path behind a default-on `parallel` feature, bit-identical to serial (test).
- [ ] Test: a uniform grid with no diffusion over a Phase-1 config gives every voxel equal to the `UnitCell` result within rtol 1e-6.

**Files to Touch:**
- `src/voxel/reaction.rs`
- `src/solver.rs`
- `Cargo.toml`

**Verification Step:**

Run:
```bash
cargo test --release voxel::reaction && cargo test --release --no-default-features voxel::reaction
```

Expected:
Tests pass with and without the parallel feature.

**Notes:**

Per-thread solver/workspace reuse. Check whether rsoda can restart cheaply; if not, note it in the commit log.
