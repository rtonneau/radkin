# Ticket 01: geant4-irt-syn-project

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `${DEV_DIR}/GEANT4/SIM/radkin` is a copy of geant4-rs (without build/ and out/) that builds against `G4_ROOT`
- [ ] The water cube is 150 µm, TimeStepModel IRT_syn, 100 events, and physics and random seeds are as in geant4-rs
- [ ] Seed CSVs (radkin format, nm) and edep per event are written at 10 ns and at 50 ns, from the same events
- [ ] The molecular reaction table and molecule definitions (name, D, charge, reaction type, k) are dumped to a text file at init

**Files to Touch:**
- `${DEV_DIR}/GEANT4/SIM/radkin/src/*`
- `${DEV_DIR}/GEANT4/SIM/radkin/macros/*`
- `${DEV_DIR}/GEANT4/SIM/radkin/CMakeLists.txt`

**Verification Step:**

Run:
```bash
pwsh ${DEV_DIR}/GEANT4/SIM/radkin/run.ps1 macros/seeds_10ns.mac ; pwsh ${DEV_DIR}/GEANT4/SIM/radkin/run.ps1 macros/seeds_50ns.mac
```

Expected:
Exit 0. out/ holds seed CSVs and edep for 100 events at each time, and a reaction-table dump. Molecule counts per event are lower at 50 ns than at 10 ns.

**Notes:**

Read the IRT_syn and scheduler sources before changing the code. Keep the edep pattern of geant4-rs.
