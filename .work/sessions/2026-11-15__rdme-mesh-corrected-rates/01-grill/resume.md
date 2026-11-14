# Session: rdme-mesh-corrected-rates

**Date:** 2026-11-14T23:22:41.194Z
**Status:** Grill phase complete

## Problem Statement

The stochastic voxel RDME (`src/ssa/rdme.rs`, `src/ssa/model.rs`) uses the macroscopic bimolecular rate constant unchanged in every voxel. As the spacing approaches the critical spacing `h*`, this underestimates encounters (Hellander, Hellander and Petzold 2012; Hellander & Hellander 2020 eq. 6-8). At the 6 nm voxels of `data/voxel_example.toml`, with `h*` = 1-4 nm for the radical recombinations, the error is a factor of 1.2-3, so the small-voxel G-values in `docs/ssa.md` are biased.

## Context & Constraints

- **Current behavior:** Each coarsening `Level` of `RdmeSetup` builds its own `SsaModel` for its voxel volume, with `c = k * prod(fixed) * (N_A V)^(1-m)`. `RdmeSetup::validity_bound` computes `h* = 0.2527 k / D` (D summed with multiplicity) and `run_rdme` warns when the fine spacing is below it. No propensity depends on `h` beyond the volume.
- **Pain point:** Bimolecular reactions are too slow on nm grids. In addition, `validity_bound` uses `k` instead of the pair constant `2k` for `2 A` reactions, so their `h*` is half the correct value.
- **Dependencies:** Keep the exact SSA, tau-leaping (`rdme_tau.rs`) and coarsening (`rdme_coarsen.rs`) unchanged in structure: they read `c` from the level's model. The LSODA, voxel deterministic and unit-cell SSA paths are unaffected. Golden tolerances stay as they are.
- **Tech stack:** Rust 2024, `src/ssa/` (model, rdme, rdme_coarsen, rdme_tau, config, voxel_sim), `toml`, `thiserror` (`RadkinError`), `log` warnings.

## Success Metrics

- A one-A-one-B association test (`A + B -> C`, uniform random start, cubic grid): with `rate_correction = "hellander"`, the mean association time matches `V/k` within its Monte Carlo error; with `"none"` it is measurably longer, consistent with Hellander 2012 eq. 1.
- Unit tests cover: the factor `1/(1 - h*/h)` on each coarsening level with its own `h`; reactions with `h <= h*` left uncorrected (warned); the effective `k*[fixed]` for fixed reactants; `kappa = 2k` for `2 A` reactions, in both the factor and `validity_bound`; config parsing and the errors (explicit `"hellander"` on an anisotropic grid or a unit-cell run).
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` are clean. `docs/ssa.md` documents the key and has the small-voxel table rerun with and without the correction.

## Architecture & Approach

New key `[ssa] rate_correction = "hellander" | "none"`, default `"hellander"` (ADR 0004). When enabled, for every mass-action reaction with two dynamic reactant slots and summed diffusion `D > 0`, on every coarsening level with cubic spacing `h`:

- `kappa` is the effective pair constant the RDME sees: `c * N_A V` converted to m^3/s, i.e. `k * prod(fixed concentrations)`, times 2 for a `2 A` reaction (its propensity `(k/N_A V) n(n-1)` is `2k/V` per pair; `D = 2 D_A`).
- `h* = 0.2527 kappa / D` (`C_3/6`, C_3 = 1.5164). This follows from eq. 6 with `1/k = 1/k_r + 1/(4 pi sigma D)`, so no reaction radius is needed.
- If `h > h*`: `c <- c / (1 - h*/h)`. Otherwise the reaction stays uncorrected on that level, with one warning per reaction (h, h*, reaction string). Coarser levels with `h > h*` are still corrected.

The factor is applied once while each `Level`'s `SsaModel` is built in `RdmeSetup::new`, so the exact SSA, tau-leaping and coarsening use it at no runtime cost. `validity_bound` uses the same `kappa` (fixes the `2 A` factor of 2). Scope rules:

- Anisotropic grid: an unset key warns and runs uncorrected; an explicit `"hellander"` is a config error.
- Unit-cell SSA: an unset key is ignored; an explicit `"hellander"` is a config error.
- Michaelis reactions and reactions with no diffusing reactant are never corrected.

## Assumptions & Trade-offs

- Assumes the macroscopic `k` of the chemistry is the Smoluchowski-Collins-Kimball constant, `1/k = 1/k_r + 1/(4 pi sigma D)`, and that `C_3` (derived for periodic cubes) holds for radkin's zero-flux grids when `L >> h`.
- Default on: every stochastic voxel result changes; `"none"` reproduces the earlier model (ADR 0004).
- Pseudo-bimolecular reactions with a fixed reactant are treated like any two-body reaction through `k*[fixed]`. They usually fall in `h <= h*` (`e_aq + O_minus + H2O`: `h*` = 75 nm) and stay uncorrected.
- Not done: non-local (neighbour-voxel) reactions below `h*`, Fange or Erban-Chapman alternatives, the hierarchical hRDME, corrections for unit-cell runs or for the deterministic voxel path.

## Open Questions

None.

## Notes

- Sources: `docs/literature/md_mistral/hellander2012/hellander2012_lite.md` (eq. 1, Theorem 1: N_steps = 1.5164 N) and `docs/literature/md_mistral/Hierarchical algorithm for the reaction-diffusion master equation_2020/..._lite.md` (eq. 6-8).
- Glossary: Critical spacing, Mesh correction. ADR 0004: correction on by default.
- From the scout seed (architecture review 2026-11-14).
