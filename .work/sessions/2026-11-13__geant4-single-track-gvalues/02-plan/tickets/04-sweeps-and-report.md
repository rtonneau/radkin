# Ticket 04: sweeps-and-report

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Whole track (event 0, bounding box + 2 µm margin) at 2 µm / 1 µm / 500 nm, and the segment (region + 1 µm margin) at 80 / 40 / 20 nm for about 5 events averaged, in both chemistries. Each dropped point has its reason recorded.
- [ ] Plots of G(t) for the main species (e_aq, OH_r, H_r, H2O2, H2, H3O_plus), 1 ps–10 µs on a log time axis: radkin at each voxel size against Geant4 IRT.
- [ ] A short report (`report.md`) with the setup, the initial voxel-size guess and how it held up, the voxel-size dependence, the agreement with IRT, the limitations, and the implications for the dose-rate follow-up. The README mean-field caveat is checked against the results.

**Files to Touch:**
- `.scratch/tests/2026-11-13__geant4-single-track-gvalues/runs/`
- `.scratch/tests/2026-11-13__geant4-single-track-gvalues/plots/`
- `.scratch/tests/2026-11-13__geant4-single-track-gvalues/report.md`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-13__geant4-single-track-gvalues/tools/plot_gt.py
```

Expected:
The PNG plots and the summary table are written. The report quotes G(t) at 1 ps, 1 ns, 100 ns, 1 µs and 10 µs for each case.

**Notes:**

Long runs go in the background. If the results show that the README's mean-field caveat has the wrong sign, correct it in ticket 1's README change as a follow-up commit.
