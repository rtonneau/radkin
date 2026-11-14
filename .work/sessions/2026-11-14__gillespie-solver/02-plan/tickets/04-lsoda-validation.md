# Ticket 04: lsoda-validation

**Acceptance Criteria:**
- [ ] A fixed-seed integration test compares the SSA ensemble mean with LSODA at the same output times on three cases: reaction-only (radical initial concentrations, self-reactions dominant), continuous beam, pulsed beam
- [ ] Criterion per species and time: `|mean − lsoda| ≤ 4·SEM_eff + bias`, where `SEM_eff` is the larger of the sample SEM and the Poisson SEM of the LSODA value, and `bias` is the documented finite-volume allowance; no golden tolerance is touched
- [ ] The test fails when the A + A factor is deliberately wrong (checked once by hand and recorded in the ticket log)
- [ ] The test runs in well under a minute in release mode
- [ ] `docs/ssa.md` records the cases, the criterion and the measured agreement
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `tests/ssa_validation.rs`
- `docs/ssa.md`

**Verification Step:**

Run:
```bash
cargo test --release --test ssa_validation
```

Expected:
3 passed, 0 failed.

**Notes:**

The LSODA `run` skips the first `T_ZERO` of an ON phase and adds its yield analytically; compare at times well after it. Choose the volume so that radical counts are in the thousands (finite-volume bias below the SEM).
