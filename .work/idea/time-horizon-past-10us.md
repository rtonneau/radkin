---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

Extend the time horizon past 10 us. O2 depletion and H2O2 yield, the usual dose-rate observables, play out over ms to s. The box homogenises in about 0.3 ms (estimate, L^2/6D), so the natural chain is IRT -> RDME -> homogeneous radkin.

- Not checked: whether a handoff from the voxel state to the unit cell exists. Time-coarsening and species levels do.
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.

## Changelog

- 2026-11-15 — capture — Extend the time horizon past 10 us (IRT -> RDME -> homogeneous) (radkin @ main)
