# Ticket 03: sweeps-and-report

**Status:** ✅ Done

## Local Test Result

- Benchmark: `run_sweep.py --only bench` (e0, g4, 10 ns, 250 nm): 6.4M voxels, 762 s, exit 0.
- `run_sweep.py --only main`: 40 runs (events 0–4 × {g4, air} × {10, 50 ns} × {250, 500 nm}), all with yields.csv and exit 0.
- `python -I tools/plot_gt.py`: 11 rows -> plots/summary.csv, and plots/{g4,air}_{10,50}ns.png.
- report.md has no pending cells. 125 nm is marked deferred.

## Review Notes

Base: 4a22fa7 (study files in .scratch, no radkin code change). Acceptance:
- Benchmark and 250/500 nm sweep: done.
- 125 nm: **deferred by the user** (approved for 10 and 50 ns, then held; finish chosen with 125 nm deferred). Not run.
- Report: G tables vs the old IRT, runtimes, the deferral and the dose-rate follow-up.

The report's numbers were checked against summary.csv. Key result: a 50 ns handover with 250 nm voxels is within about 10 % of IRT at 10 µs (e⁻aq 2.25 vs 2.02). The voxel size is not converged.

## Blockers / Challenges

- Stopping the parent pwsh of the detached sweep (so that 125 nm would not start) broke the python child's stdout. It died on its next print, after e0_air_50ns_dx250 had written its results, so its wall time is missing from index.csv. Restarted with Start-Process and direct file redirection; the existing results were skipped.
- The ticket was not registered as started at the time (the start output was truncated), so it was started at completion. Its time spent is therefore understated.
- A 125 nm run was started after the user's choice and stopped after about 20 s at the user's request. Its partial folder was deleted.

## Commits

- 4a22fa7 chore(gps): handoff 2026-11-14__geant4-radkin-irt-syn-gvalues (ticket 03 paused before 125 nm)

## Time Spent

0m (ticket-start.js to ticket-complete.js)
