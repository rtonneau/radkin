# Ticket 03: cli-resume-tests-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `radkin voxel <config> --resume <ckpt.npz> --out <dir>` works; the resumed run writes only to `--out`.
- [ ] Integration test: a run interrupted at a checkpoint and resumed matches the uninterrupted run's final state within run rtol.
- [ ] Tests for extending to a later t_end and for a mismatching grid (error).
- [ ] README documents `checkpoint_every`, `--resume`, the three uses and reading a checkpoint with `np.load`.

**Files to Touch:**
- `src/bin/radkin.rs`
- `tests/voxel_checkpoint.rs`
- `README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All tests pass, including `voxel_checkpoint`.

**Notes:**

Keep test grids small so the test runs in seconds.
