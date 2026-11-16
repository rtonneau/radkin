# Ticket 01: species-table-reaction-typing

**Effort:** high

**Acceptance Criteria:**
- [ ] `data/chem_g4.toml` and `data/chem_g4.diffusion.json` are in the repo, produced by `tools/export_g4_chem.py` from `tests/data/irt/g4_chemistry_dump.txt`; the export writes each reaction's Geant4 type as a key the config loader reads, not as a comment.
- [ ] `src/irt/species.rs` ships one table (van der Waals radius, charge, diffusion coefficient) for the water-radiolysis species under radkin's names, with values equal to the Geant4 dump; a network can override an entry.
- [ ] `src/irt/model.rs` builds, from a `SimEnv`, the pair reactions and the per-molecule reactions with: type, effective radius, reaction radius, Onsager radius, `k_dif`, `k_act`, contact probability and products.
- [ ] For the 52 pair reactions of the Geant4 network, the effective radius and the reaction radius equal the dump's values to 1e-9 relative.
- [ ] A reaction without a type gets one from `k_obs >= k_dif`; a test reports and pins how many of the 52 Geant4 reactions the rule types as Geant4 does.
- [ ] A species missing from the table makes its reactions fully diffusion-controlled with a logged warning.

**Files to Touch:**
- `src/irt/mod.rs`
- `src/irt/species.rs`
- `src/irt/model.rs`
- `src/lib.rs`
- `src/config/mod.rs`
- `src/chemistry/reactions.rs`
- `data/chem_g4.toml`
- `data/chem_g4.diffusion.json`
- `tools/export_g4_chem.py`
- `tests/data/irt/g4_chemistry_dump.txt`
- `tests/irt_model.rs`

**Verification Step:**

Run:
```bash
cargo test --release irt_model
```

Expected:
All `irt_model` tests pass, including the comparison of the 52 radii with the Geant4 dump.

**Notes:**

Geant4 formulas (11.4.3, `G4DNAMolecularReactionTable.cc`): type 0 has `Reff = k_obs / (4 pi D NA)` with `D` the sum of the two coefficients, or the single coefficient for identical reactants. Type 1 has `R = r1 + r2`; without charges `Reff = R`; with charges `Reff = -rc / (1 - exp(rc / R))` with `rc` the Onsager radius at 293.15 K and relative permittivity 80.1. `k_dif = 4 pi D Reff NA`, halved for identical reactants; `k_act = k_dif k_obs / (k_dif - k_obs)`; contact probability `Rs / (Rs + (k_dif / k_act)(Reff + Rs))` with `Rs = 0.29 nm`. Keep radkin's reaction index order and species order unchanged.
