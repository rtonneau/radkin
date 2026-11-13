# Strang splitting with implicit ADI diffusion and step-doubling control

_ADR 0002, recorded 2026-11-13 (session 2026-11-13__voxelize-diffusion)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

Phase 2 couples per-voxel chemistry with diffusion on a uniform Cartesian grid seeded by Geant4-DNA molecules. Each step runs reaction(dt/2) → diffusion(dt) → reaction(dt/2): every voxel is an independent rsoda run (rayon-parallel), diffusion is Crank-Nicolson/ADI with zero-flux boundaries, and dt is adapted by step-doubling, with steps aligned to seed-injection times. We chose this over a monolithic method-of-lines solve because it reuses the per-cell kernel unchanged, scales linearly in voxels, and avoids the stiffness limit D·dt/dx² at nm voxels.

Considered: monolithic sparse MOL (no splitting error, but memory and cost grow quickly; can be added later as a reference through ADR 0001's trait), explicit FTCS (sub-ps steps at nm voxels), FFT heat kernel (an extra dependency).
