---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

Port Geant4-DNA's IRT_syn into radkin (new `src/irt/`), starting from the 1 ps molecules, so that Geant4 is only needed up to the end of the physico-chemical stage.

- The algorithm (Geant4 11.4.3: `G4DNAIndependentReactionTimeStepper.cc`, `G4DiffusionControlledReactionModel.cc`, `G4DNAMakeReaction.cc`, about 800 lines):
  1. pair search within a cutoff radius (octree; about 500 nm from `G4IRTUtils::GetRCutOff`);
  2. reaction time per pair from the Smoluchowski survival probability: probability `Reff/r`, `t = (r - Reff)^2 / (4 D erfcinv(U/W)^2)`; partially diffusion-controlled reactions add a term built from `kact/(kact + kdif)`;
  3. pairs inside `Reff` at the start time react at once with the reaction's probability;
  4. jump to the earliest reaction, diffuse every molecule over that step, place reactants and products at sampled positions;
  5. clear the reaction set and resample every pair after each reaction step (the "synchronous" part).
- Cut at 1 ps, not at the physical stage: going earlier means porting the dissociation channels, `G4DNAWaterDissociationDisplacer`, electron-hole recombination and electron thermalisation, which buys nothing for dose rate. The 1 ps seed pass exists and matches Geant4 (G(e_aq) 4.422).
- radkin already parses the seeds as molecules with positions (`src/voxel/seed.rs`) and has the network and diffusion coefficients in `chem_g4.toml`. Not checked: whether that export carries the reaction type and effective radius IRT needs.
- What it buys:
  - inter-track reactions before the handover, which removes the ~1e7 Gy/s ceiling of independent per-track IRT;
  - an aerated spur stage: the Geant4 stepper only samples track-track pairs, not the background reactions (H2O, H3O+(B), OH-(B)) or O2 scavenging; in radkin they can be first-order exponential times;
  - one Geant4 pass instead of six, with the handover to the RDME in-process at any time;
  - a direct reference: Geant4 IRT_syn G(t) on the same seeds (statistical, not bitwise).
- Catches:
  - resampling all pairs after every reaction is heavy for hundreds of overlapping tracks; incremental resampling departs from the Geant4 algorithm;
  - IRT treats pairs as independent, which gets worse as tracks overlap densely;
  - the stepper is free-space; radkin needs its own reflective-wall rule to match the RDME box.
- Suggested order: port and validate on the 100 deaerated events; add first-order background scavenging; hand the particles over to the RDME grid. Do [rdme-reversible-mesh-correction] first, since it affects everything after the handover.
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.

## Changelog

- 2026-11-15 — capture — Port Geant4-DNA IRT_syn into radkin, from the 1 ps molecules (radkin @ main)
