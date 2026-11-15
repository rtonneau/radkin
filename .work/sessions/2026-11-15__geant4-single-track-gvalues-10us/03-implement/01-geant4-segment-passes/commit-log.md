# Ticket 01: geant4-segment-passes

**Status:** ✅ Done

## Local Test Result

`python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/g4check.py --copy`

```text
G(t) pass: 100 events, 61523.4 eV (615.2 eV/event)
edep 1ps/100ps/1ns/10ns/100ns: 100 events, identical
G [/100 eV] at the handover: seed pass vs G(t) pass
   1ps: e_aq 4.314/4.312  OH 5.370/5.367  H 1.141/1.110  H3Op 4.317/4.324  H2O2 0.046/0.050  H2 0.122/0.128
 100ps: e_aq 4.148/4.145  OH 4.849/4.845  H 0.965/0.962  H3Op 4.228/4.228  H2O2 0.215/0.215  H2 0.156/0.158
   1ns: e_aq 3.747/3.748  OH 3.989/3.987  H 0.923/0.917  H3Op 3.956/3.956  H2O2 0.465/0.463  H2 0.197/0.197
  10ns: e_aq 3.184/3.178  OH 3.329/3.326  H 0.928/0.925  H3Op 3.613/3.625  H2O2 0.584/0.587  H2 0.265/0.267
 100ns: e_aq 2.749/2.757  OH 2.937/2.939  H 0.922/0.930  H3Op 3.296/3.308  H2O2 0.621/0.629  H2 0.328/0.323
OK
```

## Review Notes

Base: bf006e8ab5c518ede631e251388a9163e90b241a (clean tree). Inline, no reviewer: every change is outside the radkin repo (the Geant4 project, which is not a git repo; backups in src.bak-20261115/ and include.bak-20261115/) or in the gitignored scratch dir, so this ticket has no repo commit.
- Geant4 `${DEV_DIR}/GEANT4/SIM/radkin`: `/rs/det/segmentXY`, `/rs/det/segmentZ` (defaults 2 / 10 um; world = twice the box). The gun is still at the -z face centre (GetSegmentZ). G(t) checkpoints are ported from geant4-rs (10 per decade, 1 ps to endTime - 1 ps, UserPreTimeStepAction), with the end-of-run counts kept. Macros `macros/seg2_{1ps,100ps,1ns,10ns,100ns,gt}.mac`. Built with nmake in the MSVC environment (vcvars64).
- 6 passes x 100 events: identical edep. The handover G agrees with the G(t) pass (all within 1 % or 2 sigma; H at 1 ps differs by 3 %, within Poisson). No (B) bulk species appear in the counts of this geometry.
- Outputs copied to `.scratch/.../g4/<pass>/`, plus chemistry_dump.txt.

## Blockers / Challenges

- The radkin Geant4 project had no G(t) checkpoint counting (only end counts); it was ported from geant4-rs.
- Mean edep is 615 eV/event (about 25 e_aq per event at 1 ps), so per-event G is noisy: 10 events pool about 6 keV.

## Commits

- bf006e8 chore(gps): plan 2026-11-15__geant4-single-track-gvalues-10us

## Time Spent

3m (ticket-start.js to ticket-complete.js)
