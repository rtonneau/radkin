# Implementation Plan

**Session:** irt-syn-port
**Date:** 2026-11-15T17:27:23.002Z
**Estimated effort:** 4-6 days

## Strategy

Build the IRT stage bottom-up in a new module `src/irt/`, one milestone of the Idea per ticket, so that every ticket leaves `cargo fmt`, `cargo clippy` and `cargo test --release` clean.

1. Data first: the species table and the typed reaction model, checked number by number against the reaction table Geant4 dumped at run time.
2. Then the sampling kernel, checked against analytic survival probabilities, with no solver around it.
3. Then the single-track solver with full resampling, which carries the acceptance test against Geant4 IRT_syn. It comes before every extension, so that the extensions are added to a solver already known to match.
4. Then the parts without a Geant4 reference, each with its limit case: background reactions, walls with image pairs and scheduled tracks, incremental resampling.
5. Then the in-process handover to the RDME, and the docs.

Facts found while planning, which shape the tickets:

- The Geant4 reference run (`${DEV_DIR}/GEANT4/SIM/radkin`, option3, IRT_syn) registers no `G4DNAScavengerProcess`: reactions 53 to 72 of its dump (with `H2O(B)`, `H3Op(B)`, `OHm(B)`) never fire. The acceptance configuration therefore runs the 52 pair reactions only, with background reactions switched off.
- `G4DNAIndependentReactionTimeStepper::Prepare` clears the reaction set at every interaction step, and `InitializeReactions` samples every pair again: this is full resampling.
- A pair found inside `Reff` reacts at once only at the start time, with the reaction's probability. Later, such a pair is recorded with the time left over from the previous pair of the same molecule (a Geant4 side effect). radkin does not copy this: such a pair gets no reaction time in that step.
- Reactants are put back at the positions they had when the pair was sampled, then moved by `G4DNAMakeReaction::UpdatePositionForReaction`; products are created at the mean of the diffusion-weighted reaction site and a uniform point on the segment between the reactants.
- The dump gives, for each reaction, the type, the effective radius and the reaction radius, and for each species the diffusion coefficient, the charge and the van der Waals radius. This is the reference of ticket 1.

## Tickets Overview

- **Ticket 1:** Species table and typed reaction model: every Geant4 reaction resolves to Geant4's type, radii and probability; the `k_obs >= k_dif` rule is reported against Geant4's types.
- **Ticket 2:** Pair sampling kernel: reaction times, contact reactions and reaction positions, checked against analytic survival probabilities.
- **Ticket 3:** Single-track solver with full resampling, `[irt]` config, `radkin irt`, G(t) outputs, and the acceptance comparison with Geant4 IRT_syn on 100 events.
- **Ticket 4:** Background reactions as exponential reaction times per molecule, with the decay and well-mixed box limit cases.
- **Ticket 5:** Reflective walls with image pairs, and tracks injected at scheduled times, with their three checks.
- **Ticket 6:** Incremental resampling, its agreement with full resampling, and the measured scaling up to about 100 tracks.
- **Ticket 7:** In-process handover from the IRT stage to the RDME in `radkin voxel`.
- **Ticket 8:** ADR 0007, `docs/irt.md`, README, AGENTS.md and tools README.

## Sequencing Rationale

- Ticket 2 needs the radii and rates of ticket 1. Ticket 3 needs both.
- Tickets 4 to 6 each add one feature to the solver of ticket 3 and are tested against it or against the unit-cell solvers. They are ordered by how much they change the event loop: background reactions add one event kind, walls change the pair search and the moves, incremental resampling replaces the event queue handling.
- Ticket 6 comes after walls and scheduled tracks, because the incremental rule must also handle image pairs and molecules that arrive later, and because the scaling measure needs many scheduled tracks.
- Ticket 7 needs scheduled tracks and walls (the IRT box is the RDME grid's box).
- Ticket 8 documents what was built and records the measured numbers.

## Risks & Mitigation

- **Risk:** G(t) does not match Geant4 within statistical error on the first run. → **Mitigation:** the match is checked on the same 1 ps molecules event by event, with the kernel already validated analytically (ticket 2); differences are traced against the Geant4 source, species by species, before anything else is built on the solver.
- **Risk:** the Geant4 reference is not free space: it ran in a 2 µm water box. → **Mitigation:** ticket 3 measures the fraction of molecules that would leave that box by 10 µs and, if it is not negligible, runs the acceptance case with the same box (walls of ticket 5 pulled forward as a plain position clamp matching Geant4's behaviour).
- **Risk:** full resampling is too slow beyond a few tracks. → **Mitigation:** the pair search uses a spatial grid from ticket 3 on; incremental resampling (ticket 6) is the working mode for many tracks.
- **Risk:** image pairs near edges and corners multiply the pair count. → **Mitigation:** images are only built for molecules within the cutoff radius of a wall; the cutoff radius shrinks with the time left to the end of the stage.
- **Risk:** one PR with eight tickets is long to review. → **Mitigation:** one commit per ticket, each with its own tests and a clean check run.

## Assumptions

- The Geant4 reference data of session `2026-11-15__geant4-single-track-gvalues-10us` (1 ps molecules, G(t) and deposited energy of 100 events) stays available under `.scratch/`; a small subset is copied into `tests/data/irt/` for the regression test.
- A species with a non-zero initial concentration is a background species in the IRT stage: its concentration is fixed, its consumption is counted and reported but not debited, and its track-born molecules are still particles.
- Zero-order sources (beam G-values) are not part of the IRT stage; a config that has both a beam and an `[irt]` table is refused.
- The acceptance tolerance is statistical: each species at each compared time within 3 standard errors of the difference, the two runs using the same 100 events but independent random streams.
