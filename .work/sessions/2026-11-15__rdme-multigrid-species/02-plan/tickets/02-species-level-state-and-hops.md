# Ticket 02: species-level-state-and-hops

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] With any species level > 0, the RDME keeps each species' counts on its own level's voxels and schedules events over the voxels of every level in use.
- [ ] Hops per species on its level with `D / h_L^2` per molecule and face (zero flux).
- [ ] Poisson fill and seed injection bin each species on its own level; `fine_field` spreads each species over the fine voxels of its own level.
- [ ] Reactions whose reactants and products all share one level fire with that level's model; a network with a cross-level reaction is refused with a clear error until ticket 03.
- [ ] With every species level 0, the existing single-level code runs and existing fixed-seed tests pass unchanged.
- [ ] Tests: diffusion-only of a level-1 species against the deterministic diffusion on the coarse grid; conservation of totals; `fine_field` units.

**Files to Touch:**
- `src/ssa/rdme.rs`
- `src/ssa/rdme_multigrid.rs`
- `src/ssa/mod.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/voxel_output.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
All checks pass, the new multi-grid tests included.

**Notes:**

Prefer a separate module (`rdme_multigrid.rs`) with its own trajectory type over branching inside `Rdme`'s hot loop, so the single-level path keeps its draw order. The species level used is `max(coarsening level, species level)`; coarsening is still off for multi-grid here (ticket 04). Record the species levels in the run metadata.
