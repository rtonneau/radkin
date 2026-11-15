# Ticket 04: analysis-report

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] G(t) is pooled over events as sum(N)/sum(E/100 eV), net of the background, with SEM over events and trajectories, for e_aq, OH, H, H2O2, H2, H3O+ (and O2⁻/HO2 aerated).
- [ ] Plots: per chemistry, the 5 handovers at 10 nm against the IRT reference; voxel convergence at 1 ps deaerated (20/10/5 nm); deaerated against aerated.
- [ ] `plots/summary.csv` and `report.md` cover the setup, results tables at 1 ns / 100 ns / 1 µs / 10 µs, a comparison with IRT and with the 11-13/11-14 deterministic results, runtimes, and recommendations for the dose-rate session.

**Files to Touch:**
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/analyse.py`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/plots/`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/analyse.py
```

Expected:
Writes the plots and summary.csv without error. At each handover time the radkin G equals the seed G.

**Notes:**

Follow the format of the 11-13 and 11-14 reports.
