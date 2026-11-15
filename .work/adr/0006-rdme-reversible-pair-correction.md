# RDME mesh correction of reversible pairs

_ADR 0006, recorded 2026-11-15 (session 2026-11-15__rdme-reversible-mesh-correction)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

Hellander mesh correction (ADR 0004) scaled only bimolecular propensities, which broke detailed balance for `<->` pairs such as water autoionisation: the pH 7 background drifted with voxel size (about 480 H3O+ and OH- lost 8, 19 and 26 over 10 us at 20, 10 and 5 nm). The yields of a voxel run also kept that drift, because they subtracted only the counts at `t_start`.

Two changes:

- **Shared factor.** Both directions of each pair declared with `<->` share the Hellander factor of the direction that is bimolecular once the fixed species (H2O) are dropped, so K_eq is kept (Fange 2010; Hellander 2012). On the multi-grid the factor is that of the level where the bimolecular direction runs. Pairs that are bimolecular on both sides keep separate factors and log the K_eq shift. Nothing is inferred from one-way reactions: a chemistry that writes `H3O_plus + OH_minus -> H2O` and `H2O -> H3O_plus + OH_minus` as two separate reactions (the Geant4 export `chem_g4.toml`) gets no shared factor and keeps the drift until its pair is declared with `<->`. Each pair with a bimolecular direction is logged at the start of a run.
- **Background run.** `[ssa] background_run` (default `true` for a voxel SSA run, `false` opts out; `true` on a unit cell is an error) runs a second ensemble with the same grid, chemistry, settings and record times, no seeds and a seed derived from `seed`. Its mean totals are subtracted in the G numerator only, so `G = (N_seeded - N_background) / (E / 100 eV)`. The `N_` columns stay gross, the standard deviations add in quadrature, and `snapshots.json` records the background seed, event counts and its noise (standard error of the mean change of every species). Snapshots are not net. The cost is a second ensemble of the same trajectory count (about 2x runtime and memory). The deterministic path already subtracts a background integration, which this brings to the RDME.

Accepted trade-offs: the background run roughly doubles runtime and memory (the opt-out exists for that); one background run adds its noise to small yields and several are not averaged; pairs are not detected among one-way reactions (a possible follow-up).
