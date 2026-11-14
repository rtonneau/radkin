# Ticket 01: voxel-bench-and-gate

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/voxel_bench.rs` runs the Geant4-seeded voxel case `--repeat N` times in-process and prints median wall time, accepted/rejected steps, voxel solves, nfev totals.
- [ ] A `--reference` mode (or tool script) runs the same case at tight control rtol (1e-6) and writes final per-voxel state; the bench compares a run against it and prints per-species domain-total relative error and max per-voxel scaled error.
- [ ] Baseline numbers (time and accuracy) recorded in `docs/perf/voxel.md`.

**Files to Touch:**
- `examples/voxel_bench.rs`
- `docs/perf/voxel.md`
- `docs/perf/README.md`

**Verification Step:**

Run:
```bash
cargo run --release --example voxel_bench -- --repeat 3
```

Expected:
Median time, step and solve counts, and accuracy vs reference printed; exit 0.

**Notes:**

Reference output goes in the session scratch dir (`.scratch/tests/2026-11-14__reduce-computational-burden`), not committed. Use the `parallel` feature as default; also report a single-thread number (`RAYON_NUM_THREADS=1`).
