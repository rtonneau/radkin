# Ticket 05: tau-leaping

**Acceptance Criteria:**
- [ ] `tau_leap = true` selects Cao–Gillespie–Petzold 2006 step selection (ε = 0.03, highest-order factors g_i) with critical reactions (fewer than 10 firings left) fired one at a time; when the leap is below 10/a0 it falls back to exact SSA steps
- [ ] Counts never go negative; a leap never crosses a beam transition or an output time
- [ ] Poisson sampler without libm calls, tested against the exact mean, variance and probabilities
- [ ] Fixed-seed test: tau-leaping and exact SSA ensemble means agree within 4 combined SEM on a large-count case, and the std columns agree within their sampling error
- [ ] The speed-up over exact SSA is measured on that case and recorded in `docs/ssa.md`, with the case where fast equilibria limit it
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/rng.rs`
- `src/ssa/tau.rs`
- `src/ssa/mod.rs`
- `src/ssa/ensemble.rs`
- `examples/ssa_bench.rs`
- `tests/ssa_tau.rs`
- `docs/ssa.md`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including `ssa_tau`.

**Notes:**

Sources are zero-order reactions and are leaped like the others. Michaelis–Menten reactions count as second order for g_i.
