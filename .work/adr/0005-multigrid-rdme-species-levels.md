# Stochastic voxel runs keep each species on its own fixed grid level

_ADR 0005, recorded 2026-11-15 (session 2026-11-15__rdme-multigrid-species)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

RDME cost is dominated by hops, and global coarsening (ADR 0003) keeps the whole grid fine as long as any species is clustered, so a numerous fast background species (aerated O2) hops on the radicals' grid. `[ssa.species_levels]` therefore puts each species on a fixed level of the existing 2^L hierarchy (Cao & Erban 2014, Erban & Winkelmann 2024): a bimolecular reaction fires in the voxel of its coarsest reactant, with the propensity (and the ADR 0004 mesh correction) of that level (Hellander & Hellander 2020, III.C.4), and a product on a finer level goes uniformly to a child voxel. Levels are set by hand, since no criterion bounds how coarse a species may go before its gradients (O2 depletion in a track) are smeared; the log only suggests `floor(log2 sqrt(D_s / D_slowest seeded))` for background species. Global coarsening stays and acts as a floor: a species sits at max(global level, own level).

Considered: automatic levels from the hRDME tolerance (it bounds the finest admissible mesh, not the coarsest); per-molecule level transfers after a diffusion time (hRDME III.B, needs per-molecule clocks that a count-based RDME does not have); making species levels and global coarsening exclusive (simpler, but loses coarsening of the seeded species once they mix).
