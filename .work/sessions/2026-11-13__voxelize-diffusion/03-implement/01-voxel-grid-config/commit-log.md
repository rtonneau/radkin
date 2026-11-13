# Ticket 01: voxel-grid-config

**Status:** ✅ Done

## Local Test Result

`cargo test --release voxel` → 7 passed; full `cargo test --release` all suites ok; clippy -D warnings clean.

## Review Notes

Base 1b557ace4a508842d77c860faf1c492d47707c31, clean tree. VoxelGrid (dims/spacing/origin, index/coords/center/locate/bounds) and VoxelState (species-major SoA, gather/scatter). Config: `[voxel.grid]` and `[voxel.diffusion]` parsed by `VoxelConfig` from the run TOML; `Config` (Phase-1 merge) untouched, so golden tests unchanged. Deviation: radkin species have no TOML entries, so the 'per-species field' is a `[voxel.diffusion]` table keyed by species label, with defaults copied from Geant4-DNA 11.4.3 molecule definitions (verified in source). Accepted inline.

## Blockers / Challenges

None

## Commits

- 83758df feat(voxel): voxel grid, SoA state and [voxel] config with Geant4-DNA diffusion defaults (ticket 01)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
