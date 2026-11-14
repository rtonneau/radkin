# Implementation Plan

**Session:** gillespie-solver
**Date:** 2026-11-14T20:06:01.197Z
**Estimated effort:** 4-5 days

## Strategy

Build the stochastic solver as a new module `src/ssa/` next to the LSODA path, which stays untouched. Start with the exact SSA on the unit cell, because it is the reference for everything else: beams, the ensemble and its outputs, then the LSODA validation that proves the propensities. Tau-leaping comes next and is checked against the exact SSA. The voxel RDME reuses the same propensity model per voxel and is built in three steps: exact next-subvolume method with seeds and outputs, then coarsening, then tau-leaping of reactions and hops.

Random numbers come from ChaCha8 streams (one stream per trajectory). The samplers (exponential, Poisson) use only IEEE basic operations and a pure-Rust logarithm, so a seed gives the same trajectory on every platform and thread count. Ensemble statistics are reduced in trajectory order, never in thread order.

## Tickets Overview

- **Ticket 1:** Unit-cell exact SSA: propensity model from SimMatrices, deterministic RNG, one seeded trajectory on a reaction-only case.
- **Ticket 2:** Beam sources: G-value × dose rate as zero-order source reactions with piecewise-constant propensities (continuous and pulsed beams).
- **Ticket 3:** Ensemble, outputs and selection: N trajectories (rayon), mean/std CSV, `solver` key, `[ssa]` table and `--solver` flag on `radkin run`.
- **Ticket 4:** LSODA validation: fixed-seed ensemble mean within k·SEM (k = 4) of LSODA on reaction-only, continuous-beam and pulsed-beam cases.
- **Ticket 5:** Opt-in tau-leaping (Cao–Gillespie–Petzold 2006) on the unit cell, checked against exact SSA, speed-up measured.
- **Ticket 6:** Voxel RDME, exact: next-subvolume method with hop events, Poisson background, Geant4 seeds as integer counts, ensemble snapshots and yields, `radkin voxel` with `solver = "ssa"`, checked against Strang-split LSODA.
- **Ticket 7:** Coarsening under RDME: Poisson dispersion test, count-summing merge, multinomial refine before re-injection.
- **Ticket 8:** Tau-leaping in the voxel RDME: reactions and hops leaped, critical channels fired exactly, checked against the exact RDME.

## Sequencing Rationale

Each ticket builds on committed work: the model and exact SSA (1) are needed by the beams (2); the ensemble (3) needs both; the validation (4) needs the ensemble and is the gate before anything approximate is added; tau-leaping (5) is checked against the exact SSA; the RDME (6) reuses the per-voxel model and the ensemble reduction; coarsening (7) and voxel tau-leaping (8) are checked against the exact fine-grid RDME of ticket 6.

## Risks & Mitigation

- **Risk:** a wrong rate-convention factor biases the LSODA comparison silently → **Mitigation:** ticket 1 tests the A + A propensity directly; ticket 4 compares the ensemble with LSODA on a self-reaction-dominated case.
- **Risk:** flaky statistical tests → **Mitigation:** fixed seeds, platform-independent samplers, k = 4 bounds with a Poisson floor on the SEM.
- **Risk:** fast equilibria make exact SSA slow → **Mitigation:** validation cases use short spans; event counts are logged; tau-leaping is opt-in and its limit on fast equilibria is documented.
- **Risk:** RDME cost on large grids → **Mitigation:** next-subvolume method with an indexed heap and incremental voxel rates; coarsening; measured event rates recorded in docs.
- **Risk:** the dispersion test coarsens too early at low counts → **Mitigation:** species totals are conserved exactly; the coarsened ensemble is compared with the fine-grid ensemble before the threshold mapping is fixed.

## Assumptions

- `rand_chacha` can be fetched from crates.io. Only `rand_chacha` (with its `rand_core` re-export) is added: the samplers are written in the crate, so the `rand` crate itself is not needed.
- Voxel runs have no beam source (as the LSODA voxel path): molecules enter by seeds only.
- The LSODA path, its outputs and the golden tests are not modified.
