---
bump: minor
floor: minor
---
- New `radkin voxel CONFIG --out DIR` command: reaction-diffusion chemistry on a uniform 3D voxel grid, seeded by Geant4-DNA molecule CSV files (`[voxel]`, `[[voxel.seed]]` and `[voxel.run]` config sections).
- Diffusion uses Crank-Nicolson ADI with zero-flux boundaries. Default diffusion coefficients come from Geant4-DNA and can be overridden per species in `[voxel.diffusion]`.
- Chemistry and diffusion are coupled by Strang splitting with adaptive step doubling. The per-voxel rsoda solves run in parallel (cargo feature `parallel`, on by default).
- Geant4-DNA molecule names are mapped through a built-in alias table, extendable in `[voxel.geant4]`. Molecules outside the grid are an error unless `voxel.clip = true`.
- Outputs: `yields.csv` (domain molecule counts, and G(t) when seeds give `energy_ev`), `.npy` concentration snapshots and a `snapshots.json` index.
- `tools/compare_geant4.py` plots G(t) against Geant4-DNA yields. `data/voxel_example.toml` is a worked example.
