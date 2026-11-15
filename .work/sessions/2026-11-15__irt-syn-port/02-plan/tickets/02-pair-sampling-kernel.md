# Ticket 02: pair-sampling-kernel

**Effort:** high

**Acceptance Criteria:**
- [ ] `src/irt/sample.rs` samples the reaction time of one pair: fully diffusion-controlled from `W = Reff / r` and `t = (r - Reff)^2 / (4 D erfcinv(U / W)^2)`; partially diffusion-controlled adds the acceptance `k_act / (k_act + k_dif)` and the second time, as Geant4 does.
- [ ] A pair inside `Reff` at the start time reacts at once with the reaction's contact probability.
- [ ] The reactant positions at the reaction time and the product positions follow `G4DNAMakeReaction`.
- [ ] `erfcinv` agrees with a reference to 1e-12 over (0, 2).
- [ ] The empirical survival probability of fully diffusion-controlled pairs matches `1 - (Reff / r) erfc((r - Reff) / sqrt(4 D t))` within statistical error at several times; the partially diffusion-controlled and ion-pair cases match their long-time reaction probability.
- [ ] The kernel does not allocate and draws from `SsaRng` only.

**Files to Touch:**
- `src/irt/sample.rs`
- `src/irt/mod.rs`
- `src/ssa/rng.rs`
- `tests/irt_sample.rs`

**Verification Step:**

Run:
```bash
cargo test --release irt_sample
```

Expected:
All `irt_sample` tests pass.

**Notes:**

`SsaRng` needs a standard normal draw; add it there with a deterministic method (no platform `libm` difference), since the SSA already avoids platform-dependent functions.
