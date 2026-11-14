# Ticket 03: fuse-half-reactions

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] In the dt/2 path of `advance_to`, the second half-reaction of the first half-step and the first half-reaction of the second are one dt/2 reaction solve (6 to 5 voxel solves per attempted step).
- [ ] Voxel tests green; accuracy gate passes.
- [ ] Bench: >= 5% faster than the previous kept state, else dropped and recorded.

**Files to Touch:**
- `src/voxel/driver.rs`
- `docs/perf/voxel.md`

**Verification Step:**

Run:
```bash
cargo test --release && cargo run --release --example voxel_bench -- --repeat 5
```

Expected:
Tests pass; bench table updated.

**Notes:**

Rannacher damping applies only to the first sub-step's diffusion; keep that. Voxel solve counters must reflect the fusion.
