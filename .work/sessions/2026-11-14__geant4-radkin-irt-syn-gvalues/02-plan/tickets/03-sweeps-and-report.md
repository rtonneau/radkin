# Ticket 03: sweeps-and-report

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Benchmark one 250 nm run, then run the 250 nm and 500 nm sweep for events 0–4 × {10, 50 ns} × {g4, air}
- [ ] 125 nm runs where the benchmark allows
- [ ] Plots and report.md with G(t) tables against the old IRT reference, runtimes, and the dose-rate follow-up

**Files to Touch:**
- `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/report.md`
- `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/plots/*`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/tools/plot_gt.py
```

Expected:
Plots and summary.csv are written, and the report has no pending cells, except for skipped 125 nm cases marked as such.

**Notes:**

Long runs go in a detached process (as before, the harness limit is 2 h).
