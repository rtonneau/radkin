# Implementation Plan

**Session:** geant4-single-track-gvalues
**Date:** 2026-11-13T23:05:34.669Z
**Estimated effort:** 1–2 days (plus overnight runs)

## Strategy

Build the cheap enablers first, then the expensive runs. (1) The radkin background-voxel change makes µm-margin grids affordable; it is the only committed code. (2) geant4-rs gets edep scoring, chemistry to 10 µs with a molecule counter, and a track segment geometry, and the two Geant4 setups are run, giving both the radkin seeds and the IRT reference. (3) The study assets are prepared in `.scratch`: the option3-matched and aerated chemistry configs, a converter from Geant4 output to radkin voxel configs, and a runtime benchmark that picks the practical sweep points. (4) The sweep runs, then plots and a report.

## Tickets Overview

- **Ticket 1:** radkin advances background voxels by copying one shared background integration per reaction half-step.
- **Ticket 2:** geant4-rs writes edep per event, runs the IRT chemistry to 10 µs with G(t) on log-spaced times, and has a track segment setup; both setups are run.
- **Ticket 3:** study configs (option3-matched deaerated, CONFIG_PAPER_01 aerated), a Geant4→radkin converter, and a runtime benchmark per voxel size.
- **Ticket 4:** voxel-size sweeps for the whole track and the segment in both chemistries, G(t) plots against Geant4 IRT, and a short report.

## Sequencing Rationale

The benchmark (ticket 3) must measure the code that the sweeps will run, so ticket 1 comes first. The converter needs the Geant4 output format with edep (ticket 2). The sweeps (ticket 4) need all of the above. Tickets 1 and 2 are independent, but they run one after the other because the checkout is shared.

## Risks & Mitigation

- **Risk:** ADI round-off breaks exact equality, so every voxel stays "active". → **Mitigation:** compare against the shared background within a tolerance far below the solver atol (relative plus absolute), and add a test where a uniform field goes through diffusion and stays background.
- **Risk:** Geant4 IRT to 10 µs on a full 100 keV track (~12k molecules) takes too long or crashes. → **Mitigation:** benchmark one event first. If it is too slow, use fewer events for the reference, or stop the reference earlier and say so in the report.
- **Risk:** the option3 reaction types (partially diffusion-controlled) differ from radkin's bulk rates. → **Mitigation:** use the observed bulk rate constants from option3's sources and document the expected differences in the report.
- **Risk:** the fine segment sweep (20 nm) needs more than one night. → **Mitigation:** the benchmark predicts runtimes; drop or crop the point, and record it in the report.

## Assumptions

- `DEV_DIR` and `G4_ROOT` are set; the Geant4 11.4 build in geant4-rs and `run.ps1` work as they do now.
- A Python environment with numpy and matplotlib is available for conversion and plots (run with `-I`).
- geant4-rs is not a git repository: its changes are verified by building and running, not committed. Study assets stay in `.scratch/tests/2026-11-13__geant4-single-track-gvalues/`.
