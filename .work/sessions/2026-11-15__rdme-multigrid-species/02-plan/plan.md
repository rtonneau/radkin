# Implementation Plan

**Session:** rdme-multigrid-species
**Date:** 2026-11-14T23:37:22.710Z
**Estimated effort:** 3-4 days

## Strategy

Build the multi-grid RDME bottom-up behind `[ssa.species_levels]`, keeping the existing single-level RDME as the path taken when no species level is set (bit-identical results). First the config and hierarchy depth, then the per-species-level state with hops and same-level reactions, then (after rebasing onto `main` with the ADR 0004 mesh correction) cross-level reactions, then global coarsening as a floor, then tau-leaping, and last the aerated-track validation and docs.

Blocker: tickets 03-06 need the per-level mesh correction (session `2026-11-15__rdme-mesh-corrected-rates`, ticket 02, unmerged). Before ticket 03 the branch is rebased onto `main` after that session merges; if it has not merged, the run stops at ticket 03.

## Tickets Overview

- **Ticket 1:** `[ssa.species_levels]` config, validation, hierarchy depth and the logged level suggestion.
- **Ticket 2:** Per-species-level RDME state: hops, Poisson fill, seed binning, `fine_field`, same-level reactions; cross-level networks refused for now.
- **Ticket 3:** Cross-level reactions (coarsest-reactant voxel, child selection, product placement), validated by the A->B and association tests (after the rebase).
- **Ticket 4:** Global coarsening as a floor for species levels, refining before seeds down to each species' level.
- **Ticket 5:** Tau-leaping with species levels.
- **Ticket 6:** Aerated track validation against the all-fine ensemble, docs and example.

## Sequencing Rationale

Config first so every later ticket is reachable from a run file. The state and hops (02) do not depend on the propensities, so they can land before the rebase. Cross-level reactions (03) need the per-level mesh-corrected `SsaModel`, hence the rebase gate. Coarsening (04) and tau-leaping (05) extend a working exact multi-grid RDME and are validated against it. The aerated track test and the docs (06) need every feature in place.

## Risks & Mitigation

- **Risk:** The mesh-correction session does not merge in time → **Mitigation:** tickets 01-02 do not need it; the run stops before 03 and resumes after the rebase.
- **Risk:** The multi-grid path changes fixed-seed results of runs without species levels → **Mitigation:** dispatch to the existing single-level code when every species level is 0; existing RDME tests must pass unchanged.
- **Risk:** Incremental rate updates across levels drift or miss a dependent voxel (a hop of a fine species changes the propensity of its coarse parent voxel) → **Mitigation:** debug-assert recomputation of the rates in tests; conservation tests over cross-level reactions, refines and injections.
- **Risk:** Tau-leaping across levels (05) is large → **Mitigation:** reuse the `rdme_tau.rs` structure with channels indexed by (level, voxel); validate against the exact multi-grid RDME with the existing criterion.
- **Risk:** The benefit is small on the bundled chemistry → **Mitigation:** ticket 06 records the measured hop reduction and wall time in `docs/ssa.md`.

## Assumptions

- `2026-11-15__rdme-mesh-corrected-rates` merges to `main` with per-level mesh-corrected `SsaModel`s in `RdmeSetup::new`.
- Diffusion coefficients come from the existing `[voxel.diffusion]` / Geant4-DNA defaults.
- The 2^L hierarchy is nested (each fine voxel lies in exactly one parent), as built by `coarse_grid`.
