# Implementation Plan

**Session:** geant4-single-track-gvalues-10us
**Date:** 2026-11-15T00:35:29.951Z
**Estimated effort:** 1-2 days (mostly compute)

## Strategy

Run the experiment in four steps. Generate Geant4 data for a new 2×2×2 µm segment, then write the radkin configs and run a pilot that sizes the compute. Run the sweep, then analyse and report. The pilot gates the sweep: it measures wall time and memory before about 12 cases × 10 events run. Artifacts live in the gitignored `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/`. Geant4 changes live in the Geant4 project, outside the repo. radkin code changes only if a bug shows up.

## Tickets Overview

- **Ticket 1:** Geant4 2×2×2 µm segment: configurable segment size, 5 IRT_syn seed passes (1 ps to 100 ns) plus the G(t) pass to 10 µs, 100 events, edep consistent across passes.
- **Ticket 2:** radkin configs (3 chemistries, RDME, grid = water box), seed conversion, and the pilot (1 event, 1 ps, 10 nm), with time, memory and the O2 level check.
- **Ticket 3:** Sweep driver and the full sweep: 10 nm × 5 handovers × 3 chemistries, plus 20 and 5 nm at 1 ps deaerated, events 0–9 × 4 trajectories.
- **Ticket 4:** Analysis: pooled G(t) with SEM, plots against IRT and the earlier deterministic results, summary CSV, and report.md with recommendations for the dose-rate session.

## Sequencing Rationale

Each ticket consumes the previous one's outputs: seeds, then configs and the pilot's sizing, then run outputs, then the report. The pilot comes before the sweep so the sweep can be resized cheaply.

## Risks & Mitigation

- **Risk:** RDME memory or time at 5 nm (6.4e7 voxels) or over 10 µs (fast acid-base equilibria fire continuously) → **Mitigation:** the pilot measures them. Limit RAYON_NUM_THREADS, then cut trajectories or events for 5 nm and record the change in the report.
- **Risk:** IRT_syn stop at 1 ps is not valid or does not match plain IRT → **Mitigation:** compare the G at the handover with the G(t) pass, as 11-14 did (within 0.5 %).
- **Risk:** a long sweep dies with its parent shell → **Mitigation:** launch it detached with `Start-Process` and file redirection, and make the driver skip cases that are already done.
- **Risk:** a radkin bug blocks the runs → **Mitigation:** fix it in the repo with a test and CI clean (AGENTS.md), on this branch.

## Assumptions

- `radkin voxel` with `[ssa]` and Geant4 seed CSVs works as in docs/ssa.md and the README; the earlier sessions' seed format still applies.
- The Geant4 build in `${DEV_DIR}/GEANT4/SIM/radkin/build` rebuilds with CMake against `G4_ROOT` 11.4.3.
- 32 GB RAM, 16 cores.
