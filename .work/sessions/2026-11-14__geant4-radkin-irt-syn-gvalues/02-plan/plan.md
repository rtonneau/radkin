# Implementation Plan

**Session:** geant4-radkin-irt-syn-gvalues
**Date:** 2026-11-14T07:29:50.282Z
**Estimated effort:** 1 day (plus overnight runs)

## Strategy

Build the Geant4 side first (project copy, IRT_syn, 10/50 ns seeds, reaction-table dump), then turn the dump into a radkin TOML. Then run a quick benchmark and the 250 nm runs, then the 500/125 nm sweep and the report. Work goes on `chore/radkin-irt-syn-gvalues`, stacked on `perf/background-voxels` (user). Study assets live in `.scratch` and the Geant4 project outside the repo. Each ticket commits its gps record, plus radkin fixes if any are needed.

## Tickets Overview

- **Ticket 1:** Geant4 project `${DEV_DIR}/GEANT4/SIM/radkin`: a copy of geant4-rs with a 150 µm box and IRT_syn, writing seeds and edep at 10 ns and 50 ns for 100 events, plus a reaction-table dump.
- **Ticket 2:** Export script that turns the dump into a radkin TOML (deaerated), and the aerated config, both loaded by radkin. Also study tooling: seed conversion, run and plot scripts with the old IRT G(t) as the reference.
- **Ticket 3:** radkin runs at 250 nm, 500 nm and 125 nm (if affordable) for events 0–4, both handovers and both chemistries, plus plots and a report.

## Sequencing Rationale

Ticket 2 needs ticket 1's dump and seeds, and ticket 3 needs ticket 2's configs and scripts.

## Risks & Mitigation

- **Risk:** IRT_syn cannot stop at 10/50 ns or does not expose molecule positions there → **Mitigation:** check G4 11.4 sources (G4DNAIRT_syn, G4Scheduler end time, chem6 / UHDR examples). The fallback is plain IRT with the same end time, reported as such.
- **Risk:** Geant4 reaction types (types 0–6, scavengers, bulk (B) species) have no direct radkin equivalent → **Mitigation:** use the observed rate constant k_obs for every type, with first-order reactions on scavengers mapped to pseudo-first-order rates. Record each mapping in the TOML comments.
- **Risk:** 125 nm whole-track grids (about 50M voxels) take too long → **Mitigation:** benchmark one 250 nm run first, and run 125 nm only for the cases that fit overnight.

## Assumptions

- The previous session's IRT G(t) files are in `.scratch/tests/2026-11-13__geant4-single-track-gvalues/g4/` and are reused as the reference.
- Geant4 11.4 is at `G4_ROOT`, and geant4-rs builds with `run.ps1`.
