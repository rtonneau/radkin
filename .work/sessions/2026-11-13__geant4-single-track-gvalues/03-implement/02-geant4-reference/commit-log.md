# Ticket 02: geant4-reference

**Status:** ✅ Done

## Local Test Result

Built with NMake in `geant4-rs/build/Release` (VS 18 vcvars64): `[100%] Built target geant4rs`.
Runs (`pwsh run.ps1 -Macro macros/<m>.mac -Out .scratch/.../g4/<m>`):
- segment_seeds: molecules 2691 rows, edep 8 rows (12 s); segment_gt: edep 8, gt 5558 rows (10 s).
- whole_seeds: molecules 1221041 rows, edep 100 (74 s); whole_gt: edep 100, gt 94462 rows (1676 s).
- edep.csv identical between the seed and G(t) passes (both setups): same events.
- Event 0 (100 keV deposited): G(OH) 5.41 at 1 ps -> 2.49 at 100 ns -> 2.15 at 10 us; G(e_aq) 4.25 -> 2.45 -> 1.95. 70/100 events deposit the full 100 keV (the others backscatter).

## Review Notes

Base: f5b1e82 (geant4-rs is not a git repo; a backup of its sources before the change is in %TEMP%/claude/g4rs-backup).
- edep per event: SteppingAction sums the edep in G4_WATER, EventAction writes edep_t<thread>.csv; run.ps1 merges molecules/edep/gt and takes -Macro/-Out (a driver.mac sets /control/macroPath so the macros can run common.mac from any output directory).
- G(t): deviation from the plan. The plan said to use G4MoleculeCounter (chem6 pattern). Its counts broke charge balance, so the counts now come from the scheduler's track lists in TimeStepAction::UserPreTimeStepAction: before each step [t, t+dt) every checkpoint in it is recorded (kills are processed at the end of ComputeTrackReaction, G4ITModelProcessor.cc:354). There are 10 checkpoints per decade, 1 ps to endTime - 1 ps. Counts are keyed by configuration user ID: the counter imbalance came from option3's bulk OHm(B)/H3Op(B) tracks (same GetName() as OHm/H3Op). The analysis excludes (B), and the track species balance charge to within 1 at 1 ps.
- Segment: /rs/det/segment true builds a 2x2x10 um G4_WATER box in a 20 um G4_Galactic world; the gun starts at the centre of its -z face. Macros: common.mac, whole_seeds/whole_gt (100 events), segment_seeds/segment_gt (8 events).

## Blockers / Challenges

G4GenericMessenger has no DeclareStaticProperty in 11.4.3 (DeclareProperty on the static works). The molecule counter was replaced as described. The build dir is build/Release (NMake), which needs the VS environment.

## Commits

- dd8051d chore: Geant4 reference runs for the single-track study (ticket 02)

## Time Spent

0m (ticket-start.js to ticket-complete.js)
