# Ticket 01: reversible-pair-factor

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] A new `ssa_*` test: a small seedless box at pH 7 (CONFIG_PAPER_01 water network) at several voxel sizes keeps the H3O_plus and OH_minus means within a few sigma of the unit-cell SSA over 10 us. It fails before the fix and passes after, on the single grid and the multi-grid.
- [ ] The declared `<->` pairs (forward and backward reaction indices) can be read from the matrices/model. Unit tests on the radiopyo network list the expected pairs.
- [ ] In `mesh_correct`, the direction that is bimolecular once fixed species are dropped sets the factor, and the other direction gets the same factor on each level, multi-grid levels included. Pairs bimolecular on both sides keep separate factors and log a line saying K_eq shifts.
- [ ] The corrected pairs are logged once at setup.
- [ ] fmt, clippy -D warnings and `cargo test --release` pass.

**Files to Touch:**
- `src/matrices.rs`
- `src/ssa/model.rs`
- `src/ssa/rdme.rs`
- `src/ssa/rdme_multigrid.rs`
- `src/ssa/voxel_sim.rs`
- `tests/ssa_rdme.rs`
- `tests/ssa_rdme_multigrid.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
All `ssa_*` tests pass, including the new seedless-background tests; clippy clean.

**Notes:**

Evidence of the bug is in `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md` and `runs/_bg_rdme/` (the 2 um box loses 8, 19 and 26 of about 480 at 20, 10 and 5 nm). Water autoionisation `H2O <-> H3O_plus + OH_minus` has H2O fixed, so its backward direction is bimolecular. Follow ADR 0004/0005/0006. The kernel must not allocate.
