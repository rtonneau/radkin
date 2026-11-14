# Ticket 01: geant4-irt-syn-project

**Status:** ✅ Done

## Local Test Result

Ran `run.ps1 -Macro macros\seeds_10ns.mac -Out out\seeds_10ns` and `... seeds_50ns.mac -Out out\seeds_50ns` in `D:\DEV\GEANT4\SIM\radkin`, as two detached processes in parallel (about 1.5–2 h). Both exited 0.
- 10 ns: molecules 726431 rows, edep 100 rows, gt 1179 rows.
- 50 ns: molecules 682897 rows, edep 100 rows, gt 1186 rows.
- `edep.csv` is identical in the two passes (same events): 46/100 events deposit the full 100 keV. In a 150 µm box many tracks escape, and some start near the lateral faces.
- Molecules in events 0–4: 10 ns [11815, 5227, 325, 6391, 6266], 50 ns [11048, 4863, 316, 6075, 5964]. Every event has fewer at 50 ns.
- `chemistry_dump.txt` is identical in both passes: 40 configurations, 72 reactions.
- Smoke (2 events, 10 ns): 230 s.

Outputs are copied to `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/g4/seeds_{10,50}ns/`.

## Review Notes

Base: 5fceb89c85b289738e92e5f47ad72fdfa4998461. No pre-existing changes in the project; the untracked files of the old session are its own.
- Copy of geant4-rs without build/ and out/. Renamed target `g4radkin`. Builds with NMake Release against `G4_ROOT` (11.4.3).
- Changes: `DetectorConstruction::GetSideLength` = 150 µm; `common.mac` sets TimeStepModel IRT_syn, with the same threads (8) and seeds (20261 113); new macros `seeds_10ns.mac`, `seeds_50ns.mac`, `smoke.mac`. Physics option4, chemistry option3 and the primary generator are unchanged.
- TimeStepAction: writes molecules (alive tracks only, water excluded) and counts at EndProcessing. The G(t) checkpoints are removed. Added `DumpChemistry` (std::call_once): species (name, user ID, D, charge, radius) and reactions (type, k_obs in both units, R, R_eff, products).
- Checked in the sources: option3's `ConstructTimeStepModel` maps IRT_syn to `G4DNAIndependentReactionTimeModel`, which steps all tracks synchronously, so the positions at endTime are valid. With a non-SBS model option3 sets type 1 on the partially diffusion-controlled reactions; the dump records that.
- No independent review: inline, outside-repo study code with direct output checks.

## Blockers / Challenges

- A heredoc escaping slip turned `\n` into raw newlines in the C++ source. Fixed by rewriting the file.
- The first chained run (10 ns then 50 ns) would have gone past the 2 h harness limit. Restarted as two detached parallel processes.
- The Geant4 project is not a git repo (like geant4-rs), so this ticket has no radkin code change.

## Commits

- 5fceb89 chore(gps): plan 2026-11-14__geant4-radkin-irt-syn-gvalues

## Time Spent

1h 08m (ticket-start.js to ticket-complete.js)
