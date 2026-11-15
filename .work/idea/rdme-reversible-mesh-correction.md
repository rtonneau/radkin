---
status: captured
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

RDME mesh correction breaks detailed balance of reversible pairs, and the RDME yields are not net of the background.

- `[ssa] rate_correction = "hellander"` scales a bimolecular propensity by 1/(1 - h*/h) but leaves its first-order reverse alone. For `H3O_plus + OH_minus -> H2O` (h* about 3.2 nm) with water autoionisation, the pH 7 background of a seedless 2 um box loses 8 / 19 / 26 of about 480 H3O+ and OH- over 10 us at 20 / 10 / 5 nm voxels. The well-mixed unit-cell SSA keeps them constant.
- Fix in the literature: scale the dissociation (reverse) propensity by the same factor so that K_eq is kept (Fange et al. 2010; Hellander 2012; CRDME 2018 eq. 30). This needs the reversible pairs identified in the network, with the fixed H2O.
- The SSA voxel `yields.csv` subtracts only the t_start counts. The deterministic path integrates the seedless background alongside (session 11-13); the RDME could run a seedless trajectory, or subtract the expected drift.
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md` ("Problem found"), runs in `runs/_bg_rdme/`.
- Related: the voxel loader keeps a dense f64 seed field (7.7 GB at 6.4e7 voxels, 15 species) next to the RDME counts; a sparse seed would halve the 5 nm memory.

## Changelog

- 2026-11-15 — capture — RDME mesh correction breaks detailed balance of reversible pairs (radkin @ test/rdme-single-track-gvalues)
