# Ticket 03: single-track-solver

**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `src/irt/` holds the particle state, a spatial-grid pair search with a cutoff radius, and the event loop with full resampling: jump to the earliest reaction, diffuse every molecule over the step, place reactants and products, resample.
- [ ] An `[irt]` table sets the resampling rule and the end or handover time; `radkin irt <config>` runs the IRT stage alone on the seeds of the config, one trajectory per random stream, and writes G(t) per species (mean and standard error) on the same record times and CSV layout as the RDME yields.
- [ ] A config with both a beam and `[irt]` is refused with a `RadkinError`.
- [ ] Molecule counts change only through reactions: a run without any reaction keeps every count, and atom balance holds on the Geant4 network.
- [ ] Acceptance: on the 100 events of the Geant4 reference, single-track, deaerated, free space, full resampling, G(t) of every species from 1 ps to 10 µs agrees with IRT_syn within 3 standard errors of the difference; `tools/compare_irt_geant4.py` produces the table and the report is saved in the session scratch directory.
- [ ] A regression test runs a small subset of the reference events with a fixed seed and pins the result.

**Files to Touch:**
- `src/irt/mod.rs`
- `src/irt/particles.rs`
- `src/irt/search.rs`
- `src/irt/solver.rs`
- `src/irt/config.rs`
- `src/irt/sim.rs`
- `src/irt/output.rs`
- `src/bin/radkin.rs`
- `src/error.rs`
- `data/irt_example.toml`
- `tests/data/irt/`
- `tests/irt_single_track.rs`
- `tools/compare_irt_geant4.py`

**Verification Step:**

Run:
```bash
cargo test --release irt_ && python -I tools/compare_irt_geant4.py --events 100
```

Expected:
All `irt_` tests pass, and the comparison prints every species within 3 standard errors at every compared time.

**Notes:**

The Geant4 reference is in `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/g4/` (`1ps/molecules.csv`, `gt/gt.csv`, `gt/edep.csv`). Its water box is 2 µm: measure how many molecules a free-space run puts outside it by 10 µs and say so in the report. Do not copy Geant4's left-over time for a pair found inside `Reff` after the start. Each event is one trajectory with its own stream; events run in parallel with rayon.
