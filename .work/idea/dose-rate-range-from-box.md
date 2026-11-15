---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

Pick the dose-rate range from the box numbers. One track deposits about 0.012 Gy in the box (615 eV in 8 um^3). By estimate (not measured), tracks start to overlap within the radical lifetime (about 0.5 um, 10 us) around 1e3-1e4 Gy/s, and 10 Gy in 1 us means about 800 tracks in the box.

- Above roughly 1e7 Gy/s, tracks begin to overlap before the 100 ns handover, so independent IRT per track stops being valid. That sets the upper limit of this approach.
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.
- Related: `multi-track-dose-rate-experiment` (uses this range), `irt-syn-port` (lifts the upper limit).

## Changelog

- 2026-11-15 — capture — Pick the dose-rate range from the box numbers (radkin @ main)
