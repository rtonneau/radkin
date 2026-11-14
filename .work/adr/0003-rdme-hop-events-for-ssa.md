# Stochastic voxel runs use an RDME with hop events, not Strang splitting

_ADR 0003, recorded 2026-11-14 (session 2026-11-14__gillespie-solver)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

The stochastic solver (`solver = "ssa"`) works on integer molecule counts, which the Strang reaction/ADI-diffusion step of ADR 0002 cannot advance. Stochastic voxel runs therefore solve the reaction-diffusion master equation: diffusion is a hop event of propensity `D/h²` per molecule and face, competing with the reactions in one event queue (next-subvolume method), with opt-in Cao-Gillespie-Petzold tau-leaping over reactions and hops. Coarsening stays available: a 2x2x2 block merges when a Poisson dispersion test passes (counts are summed), and refining splits counts multinomially.

Considered: stochastic reactions per voxel inside the existing Strang step with a multinomial diffusion step (parallel over voxels, but keeps a splitting error and a dt control that has no error estimate on noisy counts); coarsening disabled under SSA (simpler, but RDME cost on large fine grids is the main risk).
