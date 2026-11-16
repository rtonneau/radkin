# Ticket 06: incremental-resampling

**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `resampling = "incremental"` keeps the sampled reaction times and draws new ones only for pairs involving the molecules that reacted, were created or were injected; `"full"` stays the default.
- [ ] On the same seeds, G(t) of the incremental rule agrees with the full rule within statistical error for every species, single track and overlapping tracks.
- [ ] `examples/irt_scaling.rs` measures the run time against the number of overlapping tracks for both rules; the table up to about 100 tracks is saved in `docs/perf/`.
- [ ] About 100 overlapping tracks reach a 1 µs handover in minutes with the incremental rule, or the measured limit is stated.

**Files to Touch:**
- `src/irt/solver.rs`
- `src/irt/incremental.rs`
- `src/irt/config.rs`
- `examples/irt_scaling.rs`
- `docs/perf/irt-scaling.md`
- `tests/irt_incremental.rs`

**Verification Step:**

Run:
```bash
cargo test --release irt_ && cargo run --release --example irt_scaling -- --tracks 1,10,100
```

Expected:
All `irt_` tests pass, and the example prints one line per track count and rule with its wall time.

**Notes:**

With incremental resampling a molecule keeps its position until it is involved in an event; its position at a later time is drawn then, from the free diffusion since its last update, reflected at the walls.
