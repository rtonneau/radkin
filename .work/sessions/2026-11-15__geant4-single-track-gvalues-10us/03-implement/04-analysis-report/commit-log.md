# Ticket 04: analysis-report

**Status:** ✅ Done

## Local Test Result

`python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/analyse.py`

```text
17/17 cases complete
g4 1ps dx20: G(e_aq) at handover 4.422 (IRT pass at 1e-12 s: 4.422, 0-9)
g4 1ps dx10: G(e_aq) at handover 4.422 (...: 4.422)
g4 100ps dx10: 4.275 (4.275)   g4 1ns dx10: 3.870 (3.870)
g4 10ns dx10: 3.409 (3.409)    g4 100ns dx10: 3.096 (3.096)
g4 1ps dx5: 4.422 (4.422)
plots and summary.csv written to .../plots
```

## Review Notes

Base bf006e8. Inline. Every output is in the gitignored scratch dir.
- analyse.py pools G over events 0-9, with SEM from the trajectory spread, and subtracts the seedless RDME background drift of the same (chemistry, grid) (`runs/_bg_rdme/`). It writes plots (handover per chemistry, voxel convergence, chemistry comparison) and summary.csv. The IRT reference is events 0-9 and all 100.
- The plots were checked by eye: legends present, bands visible, handover curves start on the IRT seed G.
- report.md covers setup, tables at 1 ns / 100 ns / 10 us, comparison with IRT and with 11-13/11-14, aerated results, the background problem, runtimes and recommendations for the dose-rate session. All numbers were cross-checked against summary.csv.

## Blockers / Challenges

- radkin finding: the RDME pH background drifts on mesh-corrected grids (-8/-19/-26 H3O+ and OH- over 10 us at 20/10/5 nm). The mesh correction scales `H3O_plus + OH_minus` but not its first-order reverse, which breaks detailed balance. It is corrected in the analysis by seedless-run subtraction and logged as a backlog idea; no radkin change in this session (scope).
- The 10-event subset is not representative of 100 events (G(e_aq, 10 us) 2.73 against 2.14), so comparisons use events 0-9.

## Commits

- 826bd31 chore(gps): ticket 03-rdme-sweep done

## Time Spent

14m (ticket-start.js to ticket-complete.js)
