# Ticket 02: lean-voxel-solve

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] The voxel reaction step uses a dedicated integrate-to-t path: no `t_eval` storage, per-thread reuse of Workspace/sparse pattern/solver buffers, and each voxel's last internal step size carried as its next `first_step`.
- [ ] Serial and parallel paths give bit-identical results; existing voxel tests green unchanged.
- [ ] Bench: >= 5% faster than ticket-01 baseline and accuracy gate passes, else reverted and recorded as dropped in `docs/perf/voxel.md`.

**Files to Touch:**
- `src/voxel/reaction.rs`
- `src/voxel/driver.rs`
- `src/solver.rs`
- `docs/perf/voxel.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example voxel_bench -- --repeat 5
```

Expected:
All checks pass; bench shows the before/after table.

**Notes:**

Carried step sizes live in a per-voxel Vec on the driver (not persisted in checkpoints). If rsoda needs an API to reuse allocations, do it in the rsoda repo and bump `rev`. Unit-cell `solve` behaviour must not change.
