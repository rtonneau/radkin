# Ticket 02: adi-diffusion

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Douglas CN-ADI diffusion step for all mobile species with zero-flux boundaries (Thomas solver, no allocation per step).
- [ ] An optional number of backward-Euler substeps (Rannacher start).
- [ ] Test: point source vs the analytic 3D Gaussian, within a stated tolerance away from the boundaries.
- [ ] Test: total mass conserved to ~1e-12 relative; non-negativity on a delta seed with the Rannacher start.
- [ ] Test: degenerate 1D/2D grids work.

**Files to Touch:**
- `src/voxel/diffusion.rs`
- `src/voxel/mod.rs`

**Verification Step:**

Run:
```bash
cargo test --release voxel::diffusion
```

Expected:
Gaussian, mass-conservation and positivity tests pass.

**Notes:**

The Neumann boundary uses ghost-cell mirroring so the discrete operator conserves mass exactly.
