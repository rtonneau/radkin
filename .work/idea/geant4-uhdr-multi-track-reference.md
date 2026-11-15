---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

Get a Geant4 reference for the multi-track case. The UHDR example uses `G4DNAScavengerMaterial` and the Gillespie path, so it can scavenge O2 and handle several tracks. It is the only independent check available for the aerated dose-rate results.
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.
- Chemistry option3 with IRT cannot scavenge O2, so the aerated runs of PR #18 have no Geant4 reference.

## Changelog

- 2026-11-15 — capture — Get a Geant4 UHDR reference for the multi-track case (radkin @ main)
