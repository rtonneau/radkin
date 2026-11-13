# Ticket 02: geant4-reference

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] geant4-rs writes the deposited energy per event (eV) to a CSV (`event_id,edep_eV`), merged by `run.ps1`.
- [ ] A chemistry mode runs the IRT to 10 µs and writes the molecule counts per species and event at log-spaced times (about 10 per decade, 1 ps–10 µs), merged by `run.ps1`. The 1.001 ps seed molecules still come from the same events, with the same random seeds.
- [ ] A track segment setup (a macro or command) where the water fills only a 2×2×10 µm region at the track entrance in a vacuum world, the gun starts on the region's -z face, and about 5 events are run.
- [ ] Both setups (whole track: the 100 events or a subset; segment: about 5 events) have been run, and the outputs are in `.scratch/tests/2026-11-13__geant4-single-track-gvalues/g4/`.

**Files to Touch:**
- `${DEV_DIR}/GEANT4/SIM/geant4-rs/src/*.cc`, `include/*.hh` (new event/stepping actions, detector messenger, molecule counter)
- `${DEV_DIR}/GEANT4/SIM/geant4-rs/macros/*.mac`
- `${DEV_DIR}/GEANT4/SIM/geant4-rs/run.ps1`

**Verification Step:**

Run:
```bash
cmake --build ${DEV_DIR}/GEANT4/SIM/geant4-rs/build --config Release && pwsh ${DEV_DIR}/GEANT4/SIM/geant4-rs/run.ps1 <macro>
```

Expected:
The build succeeds. Each run writes molecules.csv, edep.csv (edep ≈ 100 keV for full-absorption events) and the G(t) counts. At 1 ps, G(OH) is about 4–6 molecules/100 eV, and it decreases over time.

**Notes:**

Check the Geant4 11.4 sources and the chem6 and molcounters examples for the molecule-counter API and for how to run the chemistry to 10 µs with IRT. Never guess the API. Keep the existing crash handler and the IRT TimeStepModel workaround. Use `${DEV_DIR}` and `G4_ROOT` from the environment.
