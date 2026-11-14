# Ticket 03: cross-level-reactions

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] Before any code: the branch is rebased onto `main` containing the merged mesh-correction session; if it is not merged, stop and report the ticket as blocked.
- [ ] A reaction fires in the voxel of its coarsest reactant, with finer reactants' counts summed over the children and that level's (mesh-corrected) propensity; the consumed fine molecule is picked in proportion to the child counts. Unimolecular and Michaelis reactions follow the same rule.
- [ ] A product on a finer level goes uniformly to one child voxel at its level; one on a coarser level goes to the parent voxel.
- [ ] Hops of a fine species update the propensities of the coarse reaction voxels they affect.
- [ ] Test: A on level 1 converting to B on level 0; B's mean field matches mgRDME eq. 24 within 4 SEM.
- [ ] Test: one A on level 0, one B on level 2; the mean association time matches `1 / k_meso` of the level-2 voxel within its MC error.
- [ ] Test: totals conserved through cross-level reactions; ticket 02's refusal is removed.

**Files to Touch:**
- `src/ssa/rdme_multigrid.rs`
- `src/ssa/model.rs`
- `tests/ssa_rdme_multigrid.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
All checks pass, including `ssa_rdme_multigrid` conversion and association tests.

**Notes:**

hRDME III.C.4 (eq. 10-12): the pair reacts as if both occupied the coarser voxel. mgRDME eq. 20/41: a product on a finer grid is placed uniformly among the gamma^3 children. Compute eq. 24's deterministic B field in the test with a small dense solve. A fine-species hop changes the summed count of its coarse parent voxel only when it crosses a parent boundary.
