# Ticket 01: background-voxels

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Each reaction half-step integrates the shared background vector once. A voxel whose state matches the background at the start of the half-step (within a tolerance far below the solver atol) receives that result instead of its own rsoda run.
- [ ] Voxels that are not background are integrated as before, and their results are bit-identical to the current code. The existing voxel tests pass unchanged.
- [ ] A test seeds a few voxels in a large grid and checks that the yields match a run with the skip disabled (within the solver tolerance), and that background voxels are actually skipped (a counter).
- [ ] A test runs a uniform background through diffusion and checks that it stays background.
- [ ] `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings` and `cargo test --release` are clean.

**Files to Touch:**
- `src/voxel/reaction.rs`
- `src/voxel/driver.rs`
- `src/voxel/sim.rs`
- `README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All three succeed: no fmt diff, no clippy warnings, all tests pass, including the new background-voxel tests.

**Notes:**

The background vector is the uniform initial state the grid starts from (the config's initial concentrations of the dynamic species), and it evolves over time: it is integrated over each half-step itself, which also covers aerated runs where it is not a steady state. Keep the kernel allocation-free. Mention the behaviour in the README's voxel section.
