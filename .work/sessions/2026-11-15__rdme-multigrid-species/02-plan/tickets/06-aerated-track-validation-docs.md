# Ticket 06: aerated-track-validation-docs

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] Fixed-seed test: aerated pure water, radicals seeded in the centre, O2 on level 1-2; species totals within 4 combined SEM of the all-fine ensemble.
- [ ] `docs/ssa.md` documents `[ssa.species_levels]`, the reaction and product rules, the floor with global coarsening, tau-leaping, and a table of hops, wall time and difference for O2 on levels 0, 1, 2 (100 nm voxels).
- [ ] README usage and a commented `[ssa.species_levels]` example in `data/voxel_example.toml`.
- [ ] Changelog fragment.

**Files to Touch:**
- `tests/ssa_rdme_multigrid.rs`
- `docs/ssa.md`
- `README.md`
- `data/voxel_example.toml`
- `examples/voxel_bench.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks pass; `docs/ssa.md` has the measured multi-grid table.

**Notes:**

Keep the test's trajectory count low enough for CI (a few seconds); run the larger measurement for the docs table by hand and record its setup.
