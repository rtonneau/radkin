# Ticket 06: acceptance

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] The baseline case runs with coarsening (max_coarsen_levels >= 2); species totals stay within rtol 1e-6 of the fine run and fields within the merge tolerance.
- [ ] `coarsen_tol_factor` is tuned: the default is kept or changed with the measured evidence.
- [ ] `docs/perf/voxel.md` records the times, the speedup, the level history and the keep decision per the docs/perf keep rule, and whether the local-AMR follow-up is needed.
- [ ] `README.md` documents the two new `[voxel]` keys; full CI checks are clean.

**Files to Touch:**
- `docs/perf/voxel.md`
- `examples/voxel_bench.rs`
- `README.md`
- `src/voxel/config.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example voxel_bench
```

Expected:
All checks pass; the bench prints fine and coarsened times and the totals agree within rtol 1e-6.

**Notes:**

If the speedup fails the keep rule, say so plainly in docs/perf and leave the feature off by default; don't loosen tolerances to pass.
