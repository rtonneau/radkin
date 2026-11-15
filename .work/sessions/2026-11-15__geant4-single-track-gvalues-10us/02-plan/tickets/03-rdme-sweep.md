# Ticket 03: rdme-sweep

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] The sweep driver runs every case and event, skips cases that are done, and logs wall times. It is launched detached.
- [ ] Every case is complete: 10 nm × {1 ps, 100 ps, 1 ns, 10 ns, 100 ns} × 3 chemistries, plus 20 and 5 nm at 1 ps deaerated, events 0–9. Any reduction from the pilot is recorded.

**Files to Touch:**
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/run_sweep.py`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/runs/`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/run_sweep.py --status
```

Expected:
Every case shows 10/10 events done, each with a yields.csv that reaches 10 µs.

**Notes:**

Use `Start-Process` with file redirection. Set RAYON_NUM_THREADS from the pilot's memory figure.
