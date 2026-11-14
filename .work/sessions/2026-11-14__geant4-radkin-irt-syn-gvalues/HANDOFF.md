# Handoff: geant4-radkin-irt-syn-gvalues

**Session:** 2026-11-14__geant4-radkin-irt-syn-gvalues
**Saved:** 2026-11-14T08:53:59.758Z
**Current phase:** ship
**Active ticket:** 03-sweeps-and-report

## Where I Stopped

Ticket 03: the sweep is running detached (started 2026-11-14 ~10:50 local): `run_sweep.py --only main` (20 runs at 250 nm and 20 at 500 nm, about 4–5 h), then `--only fine` (event 0, g4 chemistry, 10 and 50 ns at 125 nm, about 2 h each). Logs: `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/runs/main_sweep.log` and `fine_sweep.log`.

## Reasoning So Far

Benchmark: 250 nm on event 0 = 6.4M voxels, 762 s. 125 nm for all 20 cases (~40 h) is not affordable, so it is cut to event 0 with the g4 chemistry. Detached because of the 2 h harness limit. IRT_syn at the handover matches the old plain IRT within 0.5 %, which supports reusing the old reference.

## Next Step

Check that both logs show `exit 0` for every case, then run `python -I tools/plot_gt.py` in the scratch dir. Write report.md (G tables at 10 ns / 50 ns / 100 ns / 1 µs / 10 µs vs the old IRT, runtimes, the 125 nm cut, the dose-rate follow-up), fill the ticket 03 commit-log, `ticket-complete.js 3 --commit HEAD`, then `/gps finish` (PR into perf/background-voxels).

## Open Questions

None blocking.

## Settled Decisions (do not re-litigate)

Stacked branch chore/radkin-irt-syn-gvalues (user). 150 µm box, IRT_syn, seeds at 10 and 50 ns, Geant4-exported chemistry plus aerated, old IRT G(t) as the reference, events 0–4 pooled, previous session's sweep killed (user). The Geant4 project and study files stay outside git.

## Uncommitted Work

- **Git status (project):** clean
- **Git status (session dir):**  M .work/sessions/2026-11-14__geant4-radkin-irt-syn-gvalues/.session-config.json, ?? .work/sessions/2026-11-14__geant4-radkin-irt-syn-gvalues/03-implement/03-sweeps-and-report/commit-log.md
- **Why not committed:** Only the session record; the study files are in .scratch by design and the Geant4 project is outside the repo.

## Machine State (auto-filled)

- **Ticket queue:** 01-geant4-irt-syn-project: done, 02-chemistry-export-and-tooling: done, 03-sweeps-and-report: pending
- **Recent commits:** 5f6f3a0 chore(gps): ticket 02-chemistry-export-and-tooling done, 2112de0 chore(gps): ticket 01-geant4-irt-syn-project done, 5fceb89 chore(gps): plan 2026-11-14__geant4-radkin-irt-syn-gvalues
