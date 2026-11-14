# Ticket 02: rdme-mesh-correction

**Status:** ✅ Done

## Local Test Result

`cargo clippy --release --all-targets -- -D warnings`: clean. `cargo test --release`: all suites pass (lib 125, ssa_rdme 6, ssa_tau 3, ssa_validation 3, golden_runs 5, ...). New unit tests: `two_a_uses_the_pair_constant`, `mesh_correction_scales_each_level_with_its_own_spacing`, `mesh_correction_leaves_reactions_at_or_below_the_pole`, `mesh_correction_skips_first_order_and_still_pairs`, `mesh_correction_uses_the_effective_constant_of_a_fixed_reactant`, `mesh_correction_refuses_an_anisotropic_grid`, `correction_defaults_to_hellander_on_cubic_grids_only`.

## Review Notes

Base: 50d12577aef6df0c8cb9bb6e7a4d809dea2255b4. Pre-existing untracked `docs/literature/` is not owned.

- `critical_spacings(model, diffusion)`: `(j, h*)` for mass-action reactions of order 2 with summed `D > 0`; `kappa = c N_A V / (1e3 N_A)`, doubled for `2 A`. `h* = (1.5164/6) kappa / D` (was the rounded 0.2527; the existing test now compares at 1e-12).
- `2 A`: radiopyo stores `k = k_value / 2`; the propensity is `2c` per pair, so the pair constant is `k_value`. The old `validity_bound` used `k`, half the right `h*`. The test pins `c N_A V = 1.5e10` for `k_value = 3e10`.
- `RdmeSetup::new(..., correct)`: `h*` from the fine model (volume-independent), then every level's `c /= 1 - h*/h` where `h > h*`; anisotropic grid + correct is an error. Exposes `corrected()` and `critical_spacings()`.
- `run_rdme` / `resolve_correction`: unset → Hellander on cubic, warn + none on anisotropic; explicit Hellander on anisotropic → config error; "none" → none. With the correction on, one warning per reaction with `h <= h*` on the fine grid; the old single worst-reaction warning remains for uncorrected runs.
- Added `rate_correction` ("hellander"/"none", resolved) to `RdmeOutput` and `snapshots.json` so results record the model.
- The existing RDME integration tests pass unchanged with the default on (coarse voxels, factor near 1).

## Blockers / Challenges

The first `2 A` expectation assumed `mat.k = k_value`; a debug print showed radiopyo's `k_value / 2`, so the test now states the convention. Debug code removed.

## Commits

- 7e05026 feat: mesh-corrected bimolecular propensities in the RDME, 2A critical spacing fix (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
