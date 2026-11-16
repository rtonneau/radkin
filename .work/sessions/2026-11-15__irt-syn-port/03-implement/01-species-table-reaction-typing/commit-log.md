# Ticket 01: species-table-reaction-typing

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_model`: 7 passed, 0 failed (`tests/irt_model.rs`), plus 1 unit test in `src/irt/species.rs`.

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure).

Printed by the rule test: "rule agrees with Geant4 on 47/52 pair reactions; full in Geant4, partial by the rule: H_r + H_r, H_r + e_aq, H_r + O_atom, H_r + O_minus, e_aq + e_aq; partial in Geant4, full by the rule: none".

## Review Notes

Base: `aee437c1a80d57bb952ccff5c82e5064d80738f7`. Pre-existing changes: only the session record under `.work/` (owned by gps).

Against the Acceptance Criteria:
- `data/chem_g4.toml` and `data/chem_g4.diffusion.json` are produced by `tools/export_g4_chem.py` from `tests/data/irt/g4_chemistry_dump.txt`; each pair reaction carries `irt_type = "full" | "partial"`, read by `SimEnv::from_config` into `KReaction::irt_type` (k_reaction and arrhenius sections).
- `SpeciesTable::water()` equals the dump's 15 species (radius, charge, D); `set` overrides or adds an entry.
- `IrtModel::new` builds pair reactions (type, radii, `k_dif`, `k_act`, contact probability, products) and molecule reactions (first order, with fixed species, with a background species).
- The 52 effective radii and reaction radii equal the dump at 1e-9 relative, with CLHEP's constants for the Onsager radius.
- The rule `k_obs >= k_dif` is pinned at 47/52. The 5 differences are reactions Geant4 types fully diffusion-controlled although their rate constant is below the diffusion limit of the van der Waals radii. Geant4 never types "partial" a reaction at or above its limit.
- A species missing from the table gives a fully diffusion-controlled reaction and one warning.

Finding, fixed here: the export of session 2026-11-14 wrote Geant4's `k_obs` unchanged for the six self-reactions (`2 H_r`, `2 e_aq`, `2 OH_r`, `2 HO2`, `2 O_minus`, `O_atom + O_atom`). Geant4's `k_obs` is the rate constant of the reaction events (its `Reff` and `k_dif` for identical reactants show it), while radkin divides a scalar `k_value` of `2 A -> ...` by 2. The earlier `chem_g4.toml` therefore ran these six reactions at half Geant4's rate. The new export writes `2 k_obs` (for example `2 OH_r -> H2O2` at 1.1e10, the literature 2k). This affects the RDME and Strang comparisons of sessions 2026-11-14 and 2026-11-15 (PR #18 report).

Not needed after all: `src/config/mod.rs` (reaction tables already pass unknown keys through).

Choices: an `irt_type` of "partial" with `k_obs >= k_dif` runs fully diffusion-controlled with a warning (its `k_act` would be negative); more than two particle reactants, negative products and Michaelis-Menten reactions are refused.

## Blockers / Challenges

None.

## Commits

- 68fc2c7 feat: IRT species table and typed reaction model (ticket 01)

## Time Spent

9m (ticket-start.js to ticket-complete.js)
