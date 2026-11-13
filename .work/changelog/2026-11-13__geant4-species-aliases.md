---
bump: minor
floor: minor
---
- Geant4 seeds: every Geant4-DNA chemistry option3 species now has a built-in radkin alias (HO2, HO2-, O2, O2-, O3, O3-, O- and atomic oxygen as `O_atom`), in both name forms.
- Voxel diffusion defaults now follow Geant4-DNA option3 for all 15 water species (e.g. OH_r 2.2e-9, H2O2 2.3e-9 m^2/s); existing voxel results change.
- `data/voxel_example.toml` adds the `O_atom` species with the 10 Geant4 option3 O(3P) reactions.
