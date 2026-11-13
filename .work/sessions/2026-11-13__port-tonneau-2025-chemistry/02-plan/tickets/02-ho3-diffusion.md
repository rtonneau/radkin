# Ticket 02: ho3-diffusion

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] DEFAULT_DIFFUSION holds HO3 = 2.0e-9, with a comment that Geant4-DNA has no HO3 and the value is O3's.
- [ ] A unit test checks diffusion_of("HO3").

**Files to Touch:**
- `src/voxel/config.rs`
- `README.md`

**Verification Step:**

Run:
```bash
cargo test --release --lib voxel::config
```

Expected:
test result: ok

**Notes:**

The array length goes from 15 to 16.
