# Ticket 04: embedded-error-estimate

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] Step-doubling replaced by an embedded splitting error estimate reusing solves of the step (target about 1.5 Strang-equivalents per accepted step); step-size controller exponent matched to the estimator order.
- [ ] ADR recorded (`domain-doc.js adr`) explaining the estimator and trade-off.
- [ ] `error_decreases_with_tolerance` and other voxel tests green; accuracy gate passes.
- [ ] Bench: >= 5% faster than the previous kept state, else dropped and recorded.

**Files to Touch:**
- `src/voxel/driver.rs`
- `.work/adr/`
- `docs/perf/voxel.md`
- `README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --example voxel_bench -- --repeat 5
```

Expected:
All checks pass; final before/after table and accuracy in `docs/perf/voxel.md`.

**Notes:**

Candidate: compare Strang (R/2 D R/2) with Lie (R D, reusing the first half-reaction extended) or a symmetric-vs-reversed Strang pair; pick what passes the gate cheapest. If kept, ticket 03's fusion may become moot: remove dead code.
