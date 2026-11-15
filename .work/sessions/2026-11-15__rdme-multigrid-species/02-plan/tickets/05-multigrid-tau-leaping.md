# Ticket 05: multigrid-tau-leaping

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `tau_leap = true` works with species levels: channels are per-level hops of every species and reactions of every reaction voxel, cross-level included.
- [ ] The Cao-Gillespie-Petzold bound runs over every (species, voxel of its level) population; critical channels fire one at a time; negative counts reject and halve the leap, as today.
- [ ] Test: tau-leaping with species levels agrees with the exact multi-grid RDME within 4 combined SEM plus `tau_epsilon` times the species peak (fixed seed); standard deviations agree.

**Files to Touch:**
- `src/ssa/rdme_tau.rs`
- `src/ssa/rdme_multigrid.rs`
- `tests/ssa_rdme_multigrid.rs`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
All checks pass, including the multi-grid tau-leaping test.

**Notes:**

A cross-level reaction channel's consumption of a finer reactant is spread over the children in proportion to their counts (multinomial on the Poisson firing count). Keep the single-level tau path untouched.
