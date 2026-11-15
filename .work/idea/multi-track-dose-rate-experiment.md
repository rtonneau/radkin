---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

Run the multi-track experiment in the same 2 um box. The seed config already has a `time` and an `offset` per seed, and both solvers inject seeds at scheduled times, so no new machinery is needed for a first pass.

- Use IRT_syn seeds at a 10-100 ns handover, 10 nm voxels, and all 100 events as the track pool.
- Inject tracks at random arrival times and lateral offsets, sweeping the mean dose rate.
- Compare G(t) against the same tracks run one at a time; the difference is the inter-track effect.
- Came up on 2026-11-15 after PR #18 (RDME G(t) on single Geant4-DNA tracks), as a next step towards studying the impact of dose rate on water radiolysis (electron beam, Geant4 coupling).
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.
- Related: `dose-rate-range-from-box` (which dose rates to sweep), `rdme-reversible-mesh-correction` (prerequisite for pH-dependent species).

## Changelog

- 2026-11-15 — capture — Run the multi-track dose-rate experiment in the 2 um box (radkin @ main)
