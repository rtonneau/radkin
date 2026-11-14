# Handoff: geant4-radkin-irt-syn-gvalues

**Session:** 2026-11-14__geant4-radkin-irt-syn-gvalues
**Saved:** 2026-11-14T12:47:04.570Z
**Current phase:** ship
**Active ticket:** 03-sweeps-and-report

## Where I Stopped

Ticket 03: the 250/500 nm sweep is done (40 runs, all exit 0) and analysed (plots/summary.csv). The 125 nm runs (e0, g4 chemistry, 10 and 50 ns, ~2.5 h each) were approved, then held by the user: not started, no partial output. Work was stopped by the user, who asked to commit and push.

## Reasoning So Far

Results (pooled e0-4, deaerated, 10 µs): IRT e_aq 2.02 / OH 2.23. radkin from 10 ns: 2.69/2.82 (500 nm), 2.50/2.69 (250 nm). From 50 ns: 2.37/2.55 (500), 2.25/2.46 (250). The handover time matters most. The voxel size is not converged (500 -> 250 nm moves e_aq by -0.13; 125 nm is extrapolated to ~2.18). IRT_syn at the handover matches the old IRT within 0.5 %. Aerated: e_aq is fully scavenged by 10 µs, H2O2 is 1.06-1.13. Note: stopping the parent shell of a detached sweep kills the python child when it next prints; launch sweeps with Start-Process and direct file redirection.

## Next Step

Wait for the user to release the 125 nm runs: `Start-Process python -ArgumentList '-I','toolsun_sweep.py','--only','fine' -RedirectStandardOutput runsine_sweep.log` in the scratch dir. Then plot_gt.py, write report.md, fill the ticket 03 log, `ticket-complete.js 3 --commit HEAD`, `/gps finish` (PR chore/radkin-irt-syn-gvalues -> perf/background-voxels).

## Open Questions

When to run the 125 nm runs (approved for both 10 and 50 ns, then held by the user).

## Settled Decisions (do not re-litigate)

Stacked branch chore/radkin-irt-syn-gvalues (user). 150 µm box, IRT_syn, seeds at 10/50 ns, Geant4-exported chemistry plus aerated, old IRT reference, events 0-4 pooled. 125 nm limited to e0 with g4 chemistry. Study files in .scratch; the Geant4 project at D:/DEV/GEANT4/SIM/radkin is outside git.

## Uncommitted Work

- **Git status (project):** clean
- **Git status (session dir):** clean
- **Why not committed:** n/a

## Machine State (auto-filled)

- **Ticket queue:** 01-geant4-irt-syn-project: done, 02-chemistry-export-and-tooling: done, 03-sweeps-and-report: pending
- **Recent commits:** 5f6f3a0 chore(gps): ticket 02-chemistry-export-and-tooling done, 2112de0 chore(gps): ticket 01-geant4-irt-syn-project done, 5fceb89 chore(gps): plan 2026-11-14__geant4-radkin-irt-syn-gvalues
