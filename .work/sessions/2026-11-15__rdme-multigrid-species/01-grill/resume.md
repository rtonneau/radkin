# Session: rdme-multigrid-species

**Date:** 2026-11-14T23:30:58.436Z
**Status:** Grill phase complete

## Problem Statement

RDME cost is dominated by hops (`6 D / h^2` per molecule), and coarsening is global (`src/ssa/rdme_coarsen.rs`): one clustered species keeps the whole grid fine, or the whole grid coarsens and loses the radical structure. A numerous fast background species (aerated O2: about 150 molecules per 100 nm voxel) therefore hops on the radicals' fine grid. A multi-grid RDME (Cao & Erban 2014; Erban & Winkelmann 2024; Hellander & Hellander 2020) puts each species on its own fixed level of the existing 2^L hierarchy, so background species hop on coarse voxels while the seeded radicals keep fine resolution.

## Context & Constraints

- **Current behavior:** `RdmeSetup` builds one `Level` (grid, `SsaModel`, `1/h^2`) per coarsening level; an `Rdme` trajectory holds the counts of every species on the single current level, with one next-subvolume event queue over its voxels. Global coarsening merges every 2x2x2 block when the Poisson dispersion test passes for every species and refines to level 0 before a seed injection. Tau-leaping (`rdme_tau.rs`) builds its channels from the current level's voxels.
- **Pain point:** Background species with large counts and large `D` cost hops at the fine spacing while carrying almost no spatial structure.
- **Dependencies:** The per-level mesh correction of ADR 0004 (ticket 02 of session `2026-11-15__rdme-mesh-corrected-rates`, unmerged) must be in place before the cross-level reaction tickets: this branch is rebased onto `main` once that session merges. With no species level set, the RDME must stay bit-identical to today's (separate code path), so existing fixed-seed tests and golden tolerances do not move. LSODA, deterministic voxel and unit-cell SSA paths are unaffected.
- **Tech stack:** Rust 2024, `src/ssa/` (rdme, rdme_coarsen, rdme_tau, model, config, voxel_sim, voxel_output), `src/voxel/config.rs`, `toml`, `thiserror` (`RadkinError`), `log`, rayon over trajectories.

## Success Metrics

- Uniform first-order conversion test: A on level 1 converting to B on level 0; B's ensemble mean field matches the deterministic solution of mgRDME eq. 24 (Erban & Winkelmann 2024) within 4 SEM.
- Cross-level association test: one A on the fine level, one B on a coarser level; the mean association time matches `1 / k_meso` of the coarse voxel (hRDME eq. 12) within its Monte Carlo error.
- Aerated track test: pure aerated water, radicals seeded in the centre, O2 on level 1-2; species totals within 4 combined SEM of the all-fine ensemble (fixed seed). The hop reduction and wall time are recorded in `docs/ssa.md`.
- Tau-leaping with species levels agrees with the exact multi-grid RDME under the existing criterion (4 combined SEM plus `tau_epsilon` times the species peak).
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` clean; with no species level set, existing fixed-seed RDME tests unchanged.

## Architecture & Approach

New table `[ssa.species_levels]` (e.g. `O2 = 2`), default 0 for every species (ADR 0005). Config errors: unknown species, negative or non-integer level, use on unit-cell or deterministic runs. The hierarchy depth is `max(voxel.max_coarsen_levels, deepest species level)`, with the existing divisibility check on the grid dimensions.

- **Levels.** A species sits at `max(current coarsening level, its species level)`. Global coarsening keeps its dispersion test and cadence (global coarsening is a floor); refining before a seed stops at each species' level.
- **State.** Counts of each species on its own level's voxels; event scheduling covers the voxels of every level in use.
- **Hops.** Per species on its own level, `D / h_L^2` per molecule and face.
- **Reactions.** A reaction fires in the voxel of its coarsest reactant, with the counts of finer reactants summed over that voxel's children, and the propensity of that level's `SsaModel` (including the ADR 0004 mesh correction for that level's `h`), as in hRDME III.C.4. The consumed molecule of a finer reactant is picked in proportion to the child counts. Unimolecular and Michaelis reactions use the same coarsest-reactant rule.
- **Products.** A product on a finer level than the reaction's goes uniformly to one child voxel (mgRDME eq. 20); one on a coarser level goes to the parent voxel.
- **Initial solution and seeds.** Poisson fill and seed binning on each species' level.
- **Outputs.** `fine_field` spreads each species' count evenly over the fine voxels of its own level; species levels are recorded in the run metadata.
- **Tau-leaping.** Supported: channels are the per-level hops of every species and the reactions of every reaction voxel (cross-level included), with the Cao-Gillespie-Petzold bound over every (species, voxel-of-its-level) population.
- **Suggestion.** For species with a non-zero initial concentration, the log prints a suggested level `floor(log2 sqrt(D_s / D_slowest seeded species))`, capped by the hierarchy depth. Never applied automatically.

## Assumptions & Trade-offs

- Assumes cross-level bimolecular reactions are well described by the coarser voxel's rate (hRDME III.C.4, mean binding time argument), and that a product placed uniformly among children is acceptable (mgRDME eq. 20/25).
- Manual levels: the user owns the accuracy/cost trade-off; a too-coarse O2 level smears O2 depletion in the track. The suggestion is a hint only.
- Building directly without a measurement gate: the benefit is expected only for numerous fast background species (aerated O2 at voxels of 100 nm and up); the aerated track test records the actual hop reduction.
- Not done: automatic level choice from the hRDME tolerance `W(h) < eps`, per-molecule level transfers after a diffusion time (hRDME III.B), per-species dispersion-test coarsening, per-species ADI on the deterministic path, non-nested (overlap-weighted, eq. 25) grids.

## Open Questions

None.

## Notes

- Sources: `docs/literature/md_mistral/multigrid_RDME_2024/..._lite.md` (eq. 19-25, 41) and `docs/literature/md_mistral/hierarchical_RDME_2020/..._lite.md` (III.A-C, eq. 9-14). The literature directory is untracked in the main checkout, not in this worktree.
- From the scout seed (architecture review 2026-11-14, "Worth exploring").
- Sequencing: rebase onto `main` after `2026-11-15__rdme-mesh-corrected-rates` merges, before the reaction tickets.
- Glossary: Species level (new), Coarsening level (updated: a floor for species levels). ADR 0005.
