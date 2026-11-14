# Ticket 01: baseline

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/voxel_bench.rs` can run a representative case that homogenises (seeded, then smooth everywhere at later times) on a grid with dims divisible by 4.
- [ ] Its fine-grid wall time, voxel steps and final species totals are recorded in `docs/perf/voxel.md` as the coarsening baseline, with the exact command.

**Files to Touch:**
- `examples/voxel_bench.rs`
- `docs/perf/voxel.md`

**Verification Step:**

Run:
```bash
cargo run --release --example voxel_bench
```

Expected:
The run completes and prints timing and totals matching the numbers recorded in docs/perf/voxel.md.

**Notes:**

Read `docs/perf/README.md` for the keep rule and the existing measurement format. Reuse `data/voxel_example.toml` / `voxel_example_molecules.csv` if they homogenise; otherwise add a case. Don't change library code in this ticket.
