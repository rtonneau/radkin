# Ticket 03: radkin-reuse-solver

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `Cargo.toml` pins rsoda at the pushed `feat/lsoda-restart` commit (git `rev`, no committed `[patch]`).
- [ ] In `reaction_step_background`, each rayon worker builds one `ChemSystem` and one `Lsoda` per call and restarts it for every voxel, the background solve included; the chunked and single-pulse drivers still use `new()`.
- [ ] radkin's fmt, clippy and tests are clean; the golden results are unchanged; voxel results match the pre-change run.
- [ ] The voxel bench at 1 and 16 threads (repeated) is recorded in `docs/perf/voxel.md`; if the gain is within noise under the keep rule, the radkin code change is reverted (the rev bump may stay) and the result is recorded.

**Files to Touch:**
- `Cargo.toml`
- `Cargo.lock`
- `src/solver.rs`
- `src/voxel/reaction.rs`
- `docs/perf/voxel.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All clean; all tests pass.

**Notes:**

`mat` is fixed within one `reaction_step_background` call and the source is `DoseRateSource::None`, so a `ChemSystem<'a>` borrowing the worker's scratch can live for the whole call: no `Arc` is needed. Extend `SolveScratch` or the `run_voxels` per-worker init to hold the solver. Use `examples/voxel_bench.rs` and the method in `docs/perf/` for the bench.
