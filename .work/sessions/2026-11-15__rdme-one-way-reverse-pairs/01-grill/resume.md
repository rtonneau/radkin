# Session: rdme-one-way-reverse-pairs

**Date:** 2026-11-15T14:38:22.859Z
**Status:** Grill phase complete

## Problem Statement

The reversible-pair mesh correction (session 2026-11-15__rdme-reversible-mesh-correction, ADR 0006) only reaches pairs declared with `<->`. The Geant4 study chemistry `chem_g4.toml` writes water autoionisation as two one-way reactions, `H3O_plus + OH_minus -> H2O` and `H2O -> H3O_plus + OH_minus`, so its seedless pH 7 background still drifts (-19.4 +- 5.3 ions of about 480 over 10 us at 10 nm, against +1.3 +- 5.4 with the pair declared `<->`; report in `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/report.md`).

## Context & Constraints

- **Current behavior:** `SimMatrices::reversible_pairs` holds only the declared forward/backward pairs. `RdmeSetup` gives both directions of each pair the factor of its bimolecular direction, single grid and multi-grid. One-way reactions are each corrected on their own.
- **Pain point:** configs exported from Geant4 (and other hand-written ones) do not use `<->`, so they keep the broken detailed balance.
- **Dependencies:** ticket 01 of the previous session (`ReversiblePair`, `mesh_correct`, the multi-grid `Layout` models, `log_reversible_pairs`). Reaction index and species order must not change. Unit-cell and deterministic paths are untouched.
- **Tech stack:** Rust 2024, `src/matrices.rs`, `src/ssa/`.

## Success Metrics

- Two mass-action one-way reactions whose reactant and product sides swap exactly, once fixed species are dropped, are inferred as a reversible pair and get the shared factor.
- A reaction with more than one candidate reverse gets a warning naming them, and none of them are paired.
- One-way reactions are never paired with a direction of a declared `<->` pair (CONFIG_PAPER_01's one-way `e_aq + H2O -> H_r + OH_minus` stays unpaired). Unit tests cover the match, the ambiguous case, the mixed case and a `chem_g4.toml`-style water network. fmt, clippy and `cargo test --release` pass.

## Architecture & Approach

- At matrix build, compare the dynamic-species multisets of one-way mass-action reactions (fixed species dropped). Pair A and B when reactants(A) = products(B) and products(A) = reactants(B), and A has exactly one candidate (and B too). Add them to the reversible pairs, marked as inferred, after the declared pairs.
- The existing `RdmeSetup` classification (Shared / BothBimolecular) and the multi-grid level logic then apply unchanged.
- Always on under hellander, with no new config key. `log_reversible_pairs` says "inferred from one-way reactions" for inferred pairs, and the ambiguity warning is logged at setup.
- Amend ADR 0006 (one-way reverses are now paired) and `docs/ssa.md`, and update the glossary term "Reversible pair".

## Assumptions & Trade-offs

- Matching on dynamic species only can pair reactions written with and without H2O. The user chose this over exact sides including fixed species.
- Not done: reruns of the G-value report (unit tests only, per the user), and pairing a one-way reaction with a direction of a declared pair.

## Open Questions

None.

## Notes

This is the follow-up session on branch `fix/rdme-reversible-mesh-correction`. Both sessions go into one PR. The previous session is complete except for `/gps finish`.
