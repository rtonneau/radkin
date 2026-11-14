# Ticket 02: beam-sources

**Acceptance Criteria:**
- [ ] Each radiolytic reaction of a dynamic species is a zero-order source with propensity `yield · dose_rate(t) · N_A·V`
- [ ] The dose rate is piecewise constant between the beam transitions (`Beam::make_matrix`); at a transition the pending event is discarded and redrawn, which is exact
- [ ] Continuous and pulsed beams: the mean number of molecules produced with all other reactions switched off equals `yield · dose · N_A·V` within 4 SEM (fixed seed)
- [ ] A run without beam use (`use_beam = false`) has no source events
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/model.rs`
- `src/ssa/exact.rs`
- `tests/ssa_beam.rs`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including the continuous and pulsed source tests.

**Notes:**

`Beam::at(t)` is the dose rate; `Beam::make_matrix(stop, start)` gives `[time, dose rate]` rows. Take the dose rate of an interval at its midpoint, as `chunk_dose_rate` does.
