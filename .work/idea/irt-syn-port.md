---
status: grilled
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

radkin gets its own IRT solver (new `src/irt/`), starting from the 1 ps molecules, so that Geant4 is only needed up to the end of the physico-chemical stage.

- Why:
  - a study needs six Geant4 passes today (five seed passes plus one G(t) pass);
  - independent per-track IRT stops being valid above about 1e7 Gy/s, where tracks overlap before the handover;
  - Geant4's IRT only samples track-track pairs, so it has no aerated spur stage (no O2 scavenging, no reactions with the background H2O, H3O+, OH-);
  - the RDME cannot reproduce spur •OH chemistry from 1 ps by refining the grid (PR #18 report: •OH +40 % at 10 µs, H2O2 0.24 against 0.66).
- The Geant4 algorithm it starts from (11.4.3: `G4DNAIndependentReactionTimeStepper.cc`, `G4DiffusionControlledReactionModel.cc`, `G4DNAMakeReaction.cc`, `G4DNAMolecularReactionTable.cc`):
  1. pair search within a cutoff radius (`G4IRTUtils::GetRCutOff`);
  2. reaction time per pair from the Smoluchowski survival probability: probability `Reff/r`, `t = (r - Reff)^2 / (4 D erfcinv(U/W)^2)`; a partially diffusion-controlled reaction adds a second time built from `kact/(kact + kdif)`;
  3. a pair inside `Reff` at the start time reacts at once with the reaction's probability;
  4. jump to the earliest reaction, diffuse every molecule over that step, place reactants and products at sampled positions;
  5. clear the reaction set and resample every pair after each reaction step.
- Radii in Geant4: a fully diffusion-controlled reaction has `Reff = kobs / (4 pi D NA)`; a partially diffusion-controlled one uses the sum of the two van der Waals radii, and the Onsager radius from the charges for ion pairs.
- What radkin already has:
  - seeds parsed as molecules with positions (`src/voxel/seed.rs`);
  - the Geant4 network and diffusion coefficients (`chem_g4.toml`, `chem_g4.diffusion.json`). The export holds the reaction type as a comment only, and no radius or charge;
  - the 1 ps seed pass, which matches Geant4 (G(e_aq) 4.422).
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.
- The prerequisite `rdme-reversible-mesh-correction` is done (PRs #19 and #20, ADR 0006).

## Decisions

1. **Own solver, not a copy.** It uses Geant4's IRT sampling formulas but is designed for many overlapping tracks from the start. The Geant4 match is the acceptance test, not the design constraint.
2. **Acceptance.** In its single-track, deaerated, free-space configuration, G(t) of all species from 1 ps to 10 µs agrees with Geant4 IRT_syn on the same 100 events within statistical error.
3. **Resampling.** Two selectable rules:
   - full resampling is the reference mode and carries the Geant4 match;
   - incremental resampling is accepted when it agrees with the full mode on the same seeds.
4. **Scope.** These four parts belong to this Idea:
   - first-order background reactions (reactions with fixed-concentration species, as exponential reaction times per molecule);
   - tracks injected at scheduled times (a time and an offset per seed, as the voxel solvers accept);
   - reflective walls, consistent with the RDME box;
   - the in-process handover to the RDME.
5. **Networks.** Any radkin network runs in the IRT stage, not only Geant4 exports.
   - A reaction type given by the network (a Geant4 export) wins.
   - A reaction without a type gets one from its rate constant: near its diffusion limit it is fully diffusion-controlled, slower it is partially diffusion-controlled.
   - radkin ships one species table (radius, charge, diffusion coefficient) for the water-radiolysis species, taken from Geant4. A network may override an entry.
   - A species missing from the table makes its reactions fall back to fully diffusion-controlled, with a warning.
6. **Handover.** One fixed handover time per run, set in the config. A track scheduled to arrive after it is a config error. A long pulse therefore means a late handover or an IRT-only run.
7. **Parts without a Geant4 IRT reference** (background reactions, walls, overlapping tracks) are checked on limit cases against radkin's own solvers:
   - a first-order decay against its exponential;
   - a uniform, well-mixed box against the unit-cell SSA or LSODA result;
   - a uniform distribution that stays uniform between the walls;
   - widely spaced tracks that reproduce the single-track G(t).
   Multi-track aerated results stay provisional until the Idea `geant4-uhdr-multi-track-reference` delivers a comparison.
8. **Cut at 1 ps.** Going earlier means porting the dissociation channels, `G4DNAWaterDissociationDisplacer`, electron-hole recombination and electron thermalisation, which buys nothing for dose rate.

## Open questions

- The threshold on the rate constant that separates fully from partially diffusion-controlled reactions.
- The wall rule for a sampled reaction time: Geant4's formulas are free-space.
- How many overlapping tracks the solver must handle in reasonable time.

## Risks

- IRT treats pairs as independent. This gets worse as tracks overlap densely, and nothing in scope measures that error.
- The rate-constant rule may disagree with Geant4's own reaction types on the Geant4 network.
- Networks that do not come from Geant4 have no reference at the particle level.
- With a single handover time, long pulses force a late handover or an IRT-only run, which may be slow.
- Incremental resampling departs from the Geant4 algorithm; it is only as good as its agreement with the full mode on the cases tested.
- The scope is large for one Idea: five deliverables behind one acceptance test.

## Changelog

- 2026-11-15 — capture — Port Geant4-DNA IRT_syn into radkin, from the 1 ps molecules (radkin @ main)
- 2026-11-15 — grill
