# Ticket 01: checkpoint-format

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `Checkpoint` holds format version, t, dt, grid dims/spacing/origin, species order, state data, optional background vector, injected-seed records, config summary.
- [ ] `write(path)` writes a `.npz` (`state.npy` shaped per VoxelState layout, `background.npy`, `meta.json`) via temp file + rename; `read(path)` round-trips exactly.
- [ ] `check_compatible(grid, species)` returns `RadkinError` on mismatching dims, spacing or species order.
- [ ] `np.load` can open the file (manual check documented in notes; unit test checks the npy header format).

**Files to Touch:**
- `src/voxel/checkpoint.rs`
- `src/voxel/mod.rs`
- `Cargo.toml`
- `src/error.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release checkpoint
```

Expected:
All checkpoint tests pass, no clippy warnings.

**Notes:**

No C dependencies. VoxelRun's `dt` and `rannacher_pending` are private: add accessors/setters in driver.rs if needed.
