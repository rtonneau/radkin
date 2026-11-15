# Ticket 04: background-reactions

**Effort:** high

**Acceptance Criteria:**
- [ ] A reaction between a particle and a background species (non-zero initial concentration, or a constant species such as water) is a per-molecule event with an exponential reaction time of rate `k c`; a first-order reaction is the same with rate `k`.
- [ ] The consumption of each background species is counted and written to the outputs; it is not debited.
- [ ] `[irt]` has a switch that turns background reactions off (the acceptance configuration).
- [ ] A first-order decay follows `exp(-k t)` within statistical error.
- [ ] A uniform, well-mixed box of particles with a partially diffusion-controlled reaction matches the unit-cell SSA or LSODA result within statistical error.
- [ ] The acceptance comparison of ticket 3 is unchanged with the switch off.

**Files to Touch:**
- `src/irt/model.rs`
- `src/irt/solver.rs`
- `src/irt/config.rs`
- `src/irt/output.rs`
- `tests/irt_background.rs`

**Verification Step:**

Run:
```bash
cargo test --release irt_
```

Expected:
All `irt_` tests pass, including the decay and the well-mixed box.

**Notes:**

The well-mixed box test needs walls or a periodic image to keep the density uniform; until ticket 5 it uses a box large against the diffusion length, or is completed in ticket 5 if that is not enough. Say which in the log.
