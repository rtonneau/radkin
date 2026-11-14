# Ticket 07: rdme-coarsening

**Acceptance Criteria:**
- [ ] With `max_coarsen_levels > 0` the RDME grid drops one level when, for every species, the aggregated chi-square index of dispersion of the 2x2x2 blocks is consistent with uniformity; the threshold is set from `coarsen_tol_factor`
- [ ] Merging sums the 8 counts; hop rates and the per-voxel propensity model use the spacing and volume of the current level
- [ ] Before a seed injection the grid refines to level 0 by splitting each count multinomially over the 8 children, then may coarsen again
- [ ] Species totals are conserved exactly through coarsening and refining (test)
- [ ] Snapshots stay on the fine grid (piecewise-constant prolongation)
- [ ] Fixed-seed test: a coarsened ensemble agrees with the fine-grid ensemble on species totals within 4 combined SEM on a homogenising case; the measured event reduction and the threshold mapping are recorded in `docs/ssa.md`
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/rdme.rs`
- `src/ssa/rdme_coarsen.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/mod.rs`
- `tests/ssa_rdme.rs`
- `docs/ssa.md`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including the coarsening tests.

**Notes:**

Pearson's statistic of a block with sum S over 8 equiprobable cells has exact mean 7 and variance 14·(S−1)/S given S; blocks with S = 0 are skipped.
