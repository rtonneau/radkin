# Handoff: geant4-single-track-gvalues

**Session:** 2026-11-13__geant4-single-track-gvalues
**Saved:** 2026-11-14T01:27:32.749Z
**Current phase:** ship
**Active ticket:** 04-sweeps-and-report

## Where I Stopped

Ticket 04: every sweep point is done except the 20 nm segment runs (10 runs, 24M voxels, about 1–1.5 h each). They run as a detached process started 2026-11-14 02:22 local: `python -I tools/run_sweep.py --only fine` (log `.scratch/tests/2026-11-13__geant4-single-track-gvalues/runs/fine_sweep.log`). The report is drafted with `_pending_` cells for 20 nm.

## Reasoning So Far

Tickets 01–03 are done. The voxel-size dependence is clear without 20 nm: coarse voxels underestimate recombination, and radkin stays well above IRT at every size. The 20 nm runs are on the old radkin binary (released target is locked by the sweep). Its G columns are not background-corrected, but plot_gt.py recomputes G from the N columns minus the scaled drift of a seedless voxel (runs/_blank), validated against the fixed binary. Rejected: running the sweep as a harness background task, because of its 2 h limit.

## Next Step

Check `runs/fine_sweep.log` (all 10 `seg_e*_dx20` cases `exit 0`). Then run `python -I tools/plot_gt.py` in the scratch dir, fill the 20 nm rows of report.md (table and runtime), fill ticket 04 commit-log.md, `ticket-complete.js 4 --commit 149f00d`, and `/gps finish` (PR).

## Open Questions

None blocking. To decide for the follow-up session: seed the dose-rate runs at 1 µs (or 100 ns) from the Geant4 IRT pass instead of 1 ps (report recommendation).

## Settled Decisions (do not re-litigate)

Study assets stay in .scratch (user). Geant4 counts come from track lists keyed by user ID, not G4MoleculeCounter: bulk (B) pseudo-tracks broke charge balance. The option3 config follows the code, not the comments. G-values are net of the background chemistry (radkin fix 149f00d).

## Uncommitted Work

- **Git status (project):** clean
- **Git status (session dir):**  M .work/sessions/2026-11-13__geant4-single-track-gvalues/.session-config.json, ?? .work/sessions/2026-11-13__geant4-single-track-gvalues/03-implement/04-sweeps-and-report/commit-log.md
- **Why not committed:** Only the gps session record (.session-config.json) is uncommitted; the code is committed. Study files are in .scratch by design.

## Machine State (auto-filled)

- **Ticket queue:** 01-background-voxels: done, 02-geant4-reference: done, 03-study-setup: done, 04-sweeps-and-report: pending
- **Recent commits:** 149f00d fix: Voxel G-values net of the background chemistry (ticket 04), 8b450f4 chore(gps): ticket 03-study-setup done, 802e52b fix: O_atom self-reaction in voxel_example parsed as a species (ticket 03), ba04d41 chore(gps): ticket 02-geant4-reference done, dd8051d chore: Geant4 reference runs for the single-track study (ticket 02), f5b1e82 chore(gps): ticket 01-background-voxels done, 1d22775 perf: Share one integration for background voxels (ticket 01), 69ee0d8 chore(gps): plan 2026-11-13__geant4-single-track-gvalues
