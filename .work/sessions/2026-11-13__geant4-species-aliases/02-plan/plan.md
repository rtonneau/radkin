# Implementation Plan

**Session:** geant4-species-aliases
**Date:** 2026-11-13T22:28:15.984Z
**Estimated effort:** 3 hours

## Strategy

First make the reader know every Geant4-DNA chemistry option3 species (aliases + option3 diffusion defaults), with a unit test over the full name list. Then give voxel_example a chemistry that contains O_atom and prove, with a real one-event geant4-rs fixture, that an option3 seed loads end to end.

## Tickets Overview

- **Ticket 1:** built-in aliases for all option3 solutes and option3 diffusion defaults, with a coverage unit test.
- **Ticket 2:** O_atom species and its 10 Geant4 option3 reactions in voxel_example.toml, plus a one-event geant4-rs fixture test that seeds it.

## Sequencing Rationale

Ticket 2's fixture test needs the aliases and the O_atom diffusion default from ticket 1.

## Risks & Mitigation

- **Risk:** O + H2O rate written in the wrong H2O convention (radkin may multiply by [H2O]) → **Mitigation:** read radiopyo's handling of H2O in k_reaction (e.g. `e_aq + H2O` k_value = 19 in CONFIG_PAPER_01) before writing it; state the conversion in a TOML comment.
- **Risk:** changing DEFAULT_DIFFUSION breaks existing voxel tests that assert numbers → **Mitigation:** update only tests whose expectation derives from the defaults; never touch golden tolerances.
- **Risk:** the fixture is large → **Mitigation:** pick one event with every species, trim to the smallest such event; keep it in tests/data.

## Assumptions

- Geant4 11.4.3 sources are at ${DEV_DIR}/GEANT4/geant4-v11.4.3-source/geant4-v11.4.3; radiopyo at ${DEV_DIR}/PYTHON/radiopyo.
- geant4-rs output exists at D:\DEV\GEANT4\SIM\geant4-rs\out\molecules.csv (species,x,y,z,event_id,time_ps; nm; UTF-8).
