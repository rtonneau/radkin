# Session: geant4-species-aliases

**Date:** 2026-11-13T22:22:09.488Z
**Status:** Grill phase complete

## Problem Statement

radkin's Geant4 seed reader only knows 7 Geant4-DNA species. Output of Geant4-DNA chemistry option3 (e.g. the geant4-rs 100 keV e- run) also contains HO2, HO2-, O2, O2-, O3, O3-, O- and atomic O, so `radkin voxel` fails on it. Every option3 species must have a radkin counterpart.

## Context & Constraints

- **Current behavior:** `DEFAULT_GEANT4_ALIASES` (src/voxel/seed.rs) maps 14 names (7 species, GetName and tag forms). An unmapped name is used as a label and fails unless it is a dynamic species or in `voxel.geant4.ignore`. `DEFAULT_DIFFUSION` (src/voxel/config.rs) has the 7 option0 coefficients.
- **Pain point:** option3 seeds fail on `°O^0`, `HO_2°^0`, `O_2^0`, `O^-1`, `HO_2^-1` (and `O_2^-1`, `O_3^0`, `O_3^-1` when produced); no radkin chemistry has atomic O; diffusion defaults do not match option3.
- **Dependencies:** radiopyo configs (CONFIG_PAPER_01, advanced_water_radiolysis...) and golden tests stay untouched. Reference: Geant4 11.4.3 `G4ChemDissociationChannels_option1.cc`, `G4EmDNAChemistry_option3.cc`.
- **Tech stack:** Rust 2024, radkin voxel module, TOML configs.

## Success Metrics

- Every option3 configuration name (GetName and tag forms) either has a built-in alias or is a known non-solute (H2O^0, H2O(B), H3Op(B), OHm(B), NoneM); a unit test checks the full list.
- A one-event fixture from the geant4-rs output (all species present) seeds `data/voxel_example.toml` without error (with water ignored explicitly).
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` pass; golden tolerances untouched.

## Architecture & Approach

1. seed.rs: extend `DEFAULT_GEANT4_ALIASES` with HO_2°^0/HO2° -> HO2, HO_2^-1/HO2m -> HO2_minus, O_2^0/O2 -> O2, O_2^-1/O2m -> O2_minus, O_3^0/O3 -> O3, O_3^-1/O3m -> O3_minus, O^-1/Om -> O_minus, °O^0/Oxy -> O_atom (new label). Water and bulk (B) tags stay explicit ignores (not ignored by default). An aliased species absent from the chemistry stays a hard error.
2. config.rs: `DEFAULT_DIFFUSION` uses the option3 (option1 channels) coefficients for all 15 species (the 7 existing ones change, e.g. OH 2.2e-9, H2O2 2.3e-9, H3O+ 9.46e-9).
3. data/voxel_example.toml: add species O_atom with all 10 Geant4 option3 O(3P) reactions and their Geant4 rates; O_atom + H2O -> 2 OH_r at 1.9e3 /s pseudo-first-order (the Geant4 comment value, not the 1.0e3 /s in its code), written in radkin's H2O rate convention.
4. Tests: unit test over the full option3 name list; fixture test with one geant4-rs event.

## Assumptions & Trade-offs

- O_atom is a radkin-only label; only voxel_example gets its chemistry; other configs error on it until the user adds it or ignores it.
- Diffusion defaults now follow option3, which changes existing voxel results (expected, no golden covers them).
- Not doing: default ignore for water, selectable diffusion tables, conversion of O into 2 OH.

## Open Questions

- How radkin/radiopyo counts H2O in rate constants (CONFIG_PAPER_01 has `e_aq + H2O` k_value = 19): check radiopyo before writing the O + H2O rate.

## Notes

- Geant4 option3 O(3P) reactions: H+O->OH 2.02e10; OH+O->HO2 2.02e10; HO2+O->O2+OH 2.02e10; O+O->O2 2.20e10; H2O2+O->HO2+OH 1.60e9; H2+O->H+OH 4.77e3; OH-+O->HO2- 4.20e8; O2+O->O3 4.00e9; HO2-+O->O2-+OH 5.30e9 (M-1 s-1); O+H2O->2 OH (1.9e3 /s chosen).
- Fixture source: D:\DEV\GEANT4\SIM\geant4-rs\out\molecules.csv (option4 + chem option3, 1.001 ps, nm).
