# Ticket 02: sim-checkpoint-resume

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `[voxel.run] checkpoint_every = N` (output times, optional) and a checkpoint path/dir; a checkpoint is written every N output times and always at t_end.
- [ ] `VoxelSimulation::resume(ckpt)` starts from checkpoint t, state, dt and background, with Rannacher pending.
- [ ] Seeds already injected (per the checkpoint record) are skipped; new seeds after t_ckpt are injected; a seed before t_ckpt that was not injected is an error.
- [ ] Grid/species mismatch is a hard error; other config differences (beam, solver, t_end, seeds) are logged as warnings.
- [ ] Resuming a final checkpoint with a later t_end extends the run.

**Files to Touch:**
- `src/voxel/sim.rs`
- `src/voxel/seed.rs`
- `src/voxel/driver.rs`
- `src/voxel/checkpoint.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All tests pass.

**Notes:**

`run_schedule` rejects seeds earlier than `run.t`: filter the already-injected seeds before calling it.
