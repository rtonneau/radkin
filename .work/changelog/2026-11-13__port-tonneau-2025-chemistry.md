---
bump: minor
floor: minor
---
- New bundled config `data/tonneau_2025.toml`: the Tonneau et al. 2025 (PMB 70 235021) chemistry with the inputs of the paper's original runs (variable G-values, LET 0.56 keV/µm, a 30 Gy single pulse at 1e7 Gy/s).
- Voxel runs: HO3 now diffuses by default at 2.0e-9 m²/s. Geant4-DNA has no HO3, so it takes O3's value.
