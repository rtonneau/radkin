# Ticket 04: coarsening-floor

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] With `voxel.max_coarsen_levels > 0`, global coarsening applies to multi-grid runs: each species sits at `max(global level, species level)`.
- [ ] The dispersion test only reads species that would actually merge (species level below the next global level).
- [ ] Refining before a seed injection stops at each species' level.
- [ ] Test: conservation and `level_history` through coarsen/refine with mixed species levels; a homogenising case coarsens once mixed.

**Files to Touch:**
- `src/ssa/rdme_coarsen.rs`
- `src/ssa/rdme_multigrid.rs`
- `tests/ssa_rdme_multigrid.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
All checks pass, including the new coarsening-floor tests.

**Notes:**

Reuse the merge/split and dispersion statistics of `rdme_coarsen.rs`, restricted per species. Species already at or above the target level are untouched by a merge.
