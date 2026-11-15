---
status: grilled
scope: project
created: 2026-11-15
reviews: 0
promoted_to: null
---

## Problem

RDME mesh correction breaks detailed balance for reversible pairs. RDME yields also include the background drift.

- `[ssa] rate_correction = "hellander"` scales a bimolecular propensity by 1/(1 - h*/h) but leaves its first-order reverse unscaled. Take `H3O_plus + OH_minus -> H2O` (h* about 3.2 nm) with water autoionisation, in a seedless 2 um box at pH 7. Of about 480 H3O+ and OH-, the background loses 8, 19 and 26 over 10 us at 20, 10 and 5 nm voxels. The well-mixed unit-cell SSA keeps them constant.
- The SSA voxel `yields.csv` subtracts only the t_start counts, so any drift of the background ends up in the yields.
- Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md` ("Problem found"), runs in `runs/_bg_rdme/`.

## Decisions

- Scope: the reverse-propensity fix and background-net yields. The sparse seed field (dense f64 seed, 7.7 GB at 6.4e7 voxels) is a separate Idea.
- Reversible pairs are found automatically in the network. A bimolecular A+B -> C is matched with a first-order C -> A+B, and fixed species (H2O) are ignored. The pairs found are logged. Each matched reverse propensity is scaled by the same Hellander factor as its forward reaction, so K_eq stays the same (Fange et al. 2010; Hellander 2012; CRDME 2018 eq. 30).
- Background: a seedless RDME companion run with the same grid and settings. Its counts are subtracted from the seeded run, the same way the deterministic path integrates the seedless background. This doubles the cost.
- Acceptance: (1) a seedless box at 20/10/5 nm keeps H3O+/OH- within statistical noise of the unit-cell SSA over 10 us, as an `ssa_*` test at a small size; (2) a rerun of the single-track 10 us G-value report, compared with the deterministic voxel run.

## Open questions

- Does subtracting a seedless companion run add too much noise to small yields? Would a shared random stream or several companion runs help?
- Pairs whose product is not fixed, or that have more than one product: are they matched, and how?
- The multi-grid species-level form (ADR 0005): does it need any extra handling for the scaled reverse rate?

## Risks

- Matching pairs automatically can give a false or missed match in other networks. Logging the pairs found makes this visible.
- The companion run doubles runtime and memory, which matters at 5 nm.

## Changelog

- 2026-11-15 — capture — RDME mesh correction breaks detailed balance of reversible pairs (radkin @ test/rdme-single-track-gvalues)
- 2026-11-15 — grill
