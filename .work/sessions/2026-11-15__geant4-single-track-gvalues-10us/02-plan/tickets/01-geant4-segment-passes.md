# Ticket 01: geant4-segment-passes

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] The segment size (xy and z) can be set by macro. The default stays 2×2×10 µm, so earlier macros still work. The gun starts at the segment's −z face.
- [ ] Macros for 5 IRT_syn seed passes (endTime 1.001 ps, 100 ps, 1 ns, 10 ns, 100 ns) and one G(t) pass to 10 µs exist. All use the 2×2×2 µm segment, 100 keV e⁻, 100 events, and common.mac's threads and seeds.
- [ ] All 6 passes ran. Their outputs (molecules.csv, edep.csv, gt.csv) are copied to `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/g4/<pass>/`.
- [ ] edep per event is identical across the 6 passes. G at each handover from the seed pass agrees with the G(t) pass within 1 %.

**Files to Touch:**
- `${DEV_DIR}/GEANT4/SIM/radkin/src/DetectorConstruction.cc`
- `${DEV_DIR}/GEANT4/SIM/radkin/include/DetectorConstruction.hh`
- `${DEV_DIR}/GEANT4/SIM/radkin/src/PrimaryGeneratorAction.cc`
- `${DEV_DIR}/GEANT4/SIM/radkin/macros/seg2_*.mac`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/g4check.py`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/g4check.py
```

Expected:
edep identical across the 6 passes for all 100 events, and G per species at each handover within 1 % of the G(t) pass.

**Notes:**

Check how gt counting in the G(t) pass excludes the (B) bulk pseudo-tracks, as in 11-13. Plain IRT was used for the 11-13 1 ps seeds; here everything is IRT_syn (user decision).
