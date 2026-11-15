# Ticket 01: species-levels-config

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `[ssa.species_levels]` parses into per-species levels (default 0); unknown species, negative, non-integer or non-table values are config errors naming the key.
- [ ] Using it on a unit-cell run or a deterministic voxel run is a config error.
- [ ] Hierarchy depth is `max(voxel.max_coarsen_levels, deepest species level)`; the grid divisibility error names the axis and the key that required the depth.
- [ ] For each species with a non-zero initial concentration, the voxel SSA run logs the suggested level `floor(log2 sqrt(D_s / D_slowest seeded))`, capped by the depth; nothing is applied.
- [ ] Unit tests for parsing, errors, depth and the suggestion formula.

**Files to Touch:**
- `src/ssa/config.rs`
- `src/voxel/config.rs`
- `src/ssa/voxel_sim.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks pass; the new config tests are listed as ok.

**Notes:**

Levels are stored by dynamic species index (order of `SimMatrices`). Keep the suggestion as a pure function so it is unit-tested without a run. The table only configures; no RDME behaviour changes in this ticket.
