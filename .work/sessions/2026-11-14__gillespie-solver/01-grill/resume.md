# Session: gillespie-solver

**Date:** 2026-11-14T19:57:02.938Z
**Status:** Grill phase complete

## Problem Statement

radkin solves only the deterministic mean-field rate equations (LSODA via rsoda). At low molecule counts (small voxels, early times, low dose) the mean-field description breaks down, and there is no independent solver to cross-check rsoda. Add a second solver based on the Gillespie algorithm that takes the same inputs as the LSODA path (SimEnv reactions and species, G-values, beams, `[voxel]` config). It is validated against LSODA in the large-volume limit, then used for small-volume stochastic physics.

Promoted from the Idea `gillespie-solver` (`.work/idea/gillespie-solver.md`). Scope of this session: all six milestones of the Idea's plan.

## Context & Constraints

- **Current behavior:** `radkin` runs the unit cell and `radkin voxel` runs the voxel grid, both with LSODA. Voxel runs couple per-voxel reactions and ADI diffusion by Strang splitting (ADR 0002), with optional global coarsening (`voxel.max_coarsen_levels`). There is no RNG dependency.
- **Pain point:** no stochastic description at low counts, and no independent check of the rsoda results.
- **Dependencies:** keep the LSODA path, its outputs and the golden tests unchanged when `solver` is absent. Reuse SimEnv, the reaction and species order, G-values, beams, the `[voxel]` config, Geant4 seeds and seed aliases, and the coarsening primitives' grid pairing. radiopyo's rate convention is the reference: `2 A -> ...` has rate `k·[A]²` and stoichiometry −2.
- **Tech stack:** Rust 2024, rayon (ensemble parallelism), new dependencies `rand` and `rand_chacha` (ChaCha8Rng), `thiserror` errors (`RadkinError`), clap CLI, TOML config.

## Success Metrics

- Unit cell: with a fixed seed, the ensemble mean at large volume lies within k·SEM (k = 4) of the LSODA result on the chosen acceptance cases, for reaction-only, continuous-beam and pulsed-beam cases. The CI test is deterministic (not flaky).
- A + A reactions use the propensity `k/(N_A·V)·n·(n−1)`, removing 2 molecules per firing, and a test shows the ensemble mean matches LSODA on a second-order self-reaction.
- The same seed gives identical output regardless of thread count and platform.
- `solver = "ssa"` (config) and the CLI override select the solver; the commands are otherwise unchanged; the CSV layout matches LSODA's plus std columns.
- Tau-leaping (opt-in) matches exact SSA statistically on the unit cell and on the voxel grid, and its speed-up is measured.
- Voxel RDME: on a large-count grid, the ensemble mean agrees with the Strang-split LSODA voxel run within k·SEM; species totals are conserved exactly by hops, coarsening (sum of 8 counts) and refining (multinomial split).
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings` and `cargo test --release` stay clean; existing golden tolerances are untouched.

## Architecture & Approach

- **Selection and config:** top-level `solver = "lsoda" | "ssa"` (default `lsoda`), overridable by a CLI flag. An `[ssa]` table holds `volume` (litres; unit cell only, converts mol/L to molecule counts), `trajectories`, `seed` and `tau_leap` (default false). Voxel mode uses the voxel volume from `[voxel]`.
- **Propensities from SimEnv:** mass-action reactions use `k/(N_A·V)^(order−1)` times the falling-factorial product of reactant counts (A + A: `n·(n−1)`, matching rate `k·[A]²` with stoichiometry −2). Michaelis–Menten reactions use their rate law as the propensity. "Negative" products (the LSODA `c/(c+ε)` depletion switch) become a hard condition: the reaction fires only if enough of that species is present. Fixed species keep constant counts.
- **Beams:** each G-value × dose rate is a zero-order source reaction with a time-dependent propensity; pulses are handled by SSA with time-varying propensities.
- **Exact SSA (unit cell):** direct / next-reaction method, the reference for everything else.
- **RNG:** `rand` + `rand_chacha` ChaCha8Rng; one config seed, trajectory i uses stream i.
- **Initial counts:** unit cell: floor of `c·N_A·V` plus a Bernoulli draw of the fraction. Voxel grid: a Poisson draw per voxel with mean `c·N_A·V`. Geant4 seeds are exact integer counts binned into voxels.
- **Ensemble and output:** N trajectories in parallel with rayon, each serial; output is the ensemble mean and std on the same time grid, in the LSODA CSV layout plus std columns. Voxel output: snapshots of mean and std on the fine grid, in the existing layout.
- **Tau-leaping (opt-in):** Cao–Gillespie–Petzold 2006 step selection, critical reactions fired one at a time by exact SSA, no negative counts. Applies to the unit cell and to the voxel RDME (hops are leaped too, with negative-count guards across voxel faces).
- **Voxel RDME (ADR 0003):** diffusion is a hop event of propensity `D/h²` per molecule and face with zero-flux boundaries, competing with reactions in one event queue (next-subvolume method), in place of Strang splitting.
- **Coarsening under RDME:** supported. A 2x2x2 block passes when, for every species, its 8 counts are consistent with a uniform distribution (chi-square index of dispersion, threshold set from `coarsen_tol_factor`). Merged count = sum of the 8. Refining to level 0 before a re-injection splits each count multinomially over the 8 children. Hop rates use the spacing of the current level.
- **Validation:** fixed seed, k = 4. Unit cell against LSODA at large volume; tau-leaping against exact SSA; voxel RDME against Strang-split LSODA on a large-count grid.

Milestones (from the Idea's plan):

1. Unit-cell exact SSA: propensities from SimEnv; a seeded single trajectory on a reaction-only case.
2. Beam sources: G-value × dose rate as time-dependent zero-order propensities; pulsed and continuous beams.
3. Ensemble and output: N seeded trajectories (rayon), CSV with mean and std columns, `solver` config key and CLI flag, volume key.
4. LSODA validation: large-volume ensemble mean within k·SEM on chosen acceptance cases; fixed-seed CI test that is not flaky.
5. Opt-in tau-leaping: handles the stiff fast equilibria, matches exact SSA statistically, speed-up measured.
6. Voxel RDME: hop events at D/h², Geant4 seeds, coarsening by the dispersion test, tau-leaping of reactions and hops; checked against Strang-split LSODA on a large-count grid.

## Assumptions & Trade-offs

- Explicit CGP tau-leaping is limited by the fast equilibria (e.g. H3O+ + OH−), so its speed-up there is modest; implicit tau-leaping was rejected because it damps the fluctuations of fast species.
- Exact SSA may spend almost all its events on fast equilibria and be slow; it is still the reference.
- Each trajectory is serial; parallelism is across trajectories only.
- When the voxel spacing is below the reaction-validity bound of the RDME, a warning is logged; it is not an error.
- A wrong rate-convention factor would bias the LSODA comparison silently; the A + A test guards it.
- Statistical tests use fixed seeds and k·SEM bounds so they cannot be flaky.
- Not done: local (octree) refinement, implicit tau-leaping, any change to the LSODA path or its tolerances.

## Open Questions

- The numeric mapping from `coarsen_tol_factor` to the dispersion-test threshold, to be tuned on a homogenising case.
- The value of the reaction-validity bound on `h` used for the small-voxel warning.
- Which acceptance cases are used for the LSODA validation (chosen in the plan).

## Notes

- Glossary terms added: SSA, Propensity, Trajectory, Ensemble, Tau-leaping, Hop.
- ADR 0003 records the RDME choice over Strang splitting for stochastic voxel runs.
- Tau-leaping in the voxel RDME needs its own leap condition for hops and guards against negative counts across faces; this is the largest piece of work in milestones 5–6.
