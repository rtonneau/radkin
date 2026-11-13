# Ticket 01: option3-aliases-diffusion

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `DEFAULT_GEANT4_ALIASES` adds, in GetName and tag forms: `HO_2°^0`/`HO2°` -> HO2, `HO_2^-1`/`HO2m` -> HO2_minus, `O_2^0`/`O2` -> O2, `O_2^-1`/`O2m` -> O2_minus, `O_3^0`/`O3` -> O3, `O_3^-1`/`O3m` -> O3_minus, `O^-1`/`Om` -> O_minus, `°O^0`/`Oxy` -> O_atom; the existing 14 stay.
- [ ] Water and bulk tags are not ignored by default; an aliased species absent from the chemistry is still a hard error.
- [ ] `DEFAULT_DIFFUSION` holds the option3 (G4ChemDissociationChannels_option1 / molecule definitions as overridden there) coefficients for all 15 species, with the doc comment citing those sources.
- [ ] A unit test lists every option3 configuration name (GetName and tag forms, 17 configurations) and asserts each maps to a radkin label or is in a documented non-solute list (H2O^0, H2O(B), H3Op(B), OHm(B), NoneM).
- [ ] fmt, clippy -D warnings and the release test suite pass.

**Files to Touch:**
- `src/voxel/seed.rs`
- `src/voxel/config.rs`
- `README.md` (only if it documents the alias table or diffusion defaults)

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks pass; the new coverage test runs and passes.

**Notes:**

Verify every name and coefficient against Geant4 11.4.3 source (G4ChemDissociationChannels_option1.cc ConstructMolecule and the G4*.cc molecule definitions); GetName() is `<definition name>^<charge>` (e.g. G4O2 is "O_2", G4HO2 "HO_2°", G4Oxygen "°O", HO2m definition "HO_2" charge -1, Om definition "O" charge -1). geant4-rs output confirms: °O^0, HO_2°^0, O_2^0, O^-1, HO_2^-1. Where option1 sets a configuration-level D (H3Op 9.46e-9, °OH 2.2e-9, OHm 5.3e-9, H2 4.8e-9, H2O2 2.3e-9, HO2m 1.4e-9, ...) use that; otherwise the definition's D.
