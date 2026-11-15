# Implementation Plan

**Session:** rdme-reversible-mesh-correction
**Date:** 2026-11-15T11:58:07.935Z
**Estimated effort:** 2 days

## Strategy

Write the failing seedless-box test first. Then expose the declared `<->` pairs from the network and apply the shared Hellander factor in `mesh_correct`. Single grid and multi-grid both build on `RdmeSetup` levels, so one change there covers both, and the multi-grid level choice gets its own test. Then add the seedless background run and the net yields. Last, rerun the G-value report and write the docs.

## Tickets Overview

- **Ticket 1:** declared reversible pairs are exposed from `SimMatrices`/`SsaModel`. Both directions of each pair share the factor of their bimolecular direction in `mesh_correct` (single grid and multi-grid). Seedless-box `ssa_*` tests pass.
- **Ticket 2:** a seedless background run (default on, `[ssa] background_run = false`) makes voxel SSA `yields.csv` background-net, and the run plus its noise are recorded in the output.
- **Ticket 3:** the single-track 10 us G-value report is rerun and compared with the deterministic voxel run. `docs/ssa.md` and ADR 0006 are updated.

## Sequencing Rationale

The background run would only hide the drift if the correction were still wrong, so the correction lands first and its test shows a flat background without any subtraction. The report rerun needs both.

## Risks & Mitigation

- **Risk:** the multi-grid level picks the wrong factor for the unimolecular direction → **Mitigation:** a test on `ssa_rdme_multigrid` with species levels that checks the factor per level.
- **Risk:** a statistical test is flaky → **Mitigation:** fixed seeds, a small box, and a tolerance of several sigma from the unit-cell SSA.
- **Risk:** the background run doubles the cost of the large 5 nm runs → **Mitigation:** the opt-out flag, and the cost is noted in the docs.

## Assumptions

- The parser's forward/backward pair struct keeps the indices of both directions in the reaction order (`src/chemistry/env.rs` sets `backward.index = idx + 1`).
- The unit-cell and deterministic paths and the golden tolerances are untouched.
