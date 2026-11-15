---
bump: minor
floor: minor
---
- New `radkin irt CONFIG --out DIR` command: the chemistry on the molecules of the seeds as particles, with independent reaction times per pair (the IRT method, with Geant4-DNA's reaction-time formulas), from the molecules at 1 ps. Config section `[irt]`, with `[[irt.seed]]` as `[[voxel.seed]]`. Outputs `yields.csv` (the columns of the RDME yields) and `irt.json`. `data/irt_example.toml` is a worked example; the method and its validation are in `docs/irt.md`.
- On 100 single-track events G(t) of every species matches Geant4-DNA IRT_syn within its statistical error from 10 ps to 10 us (`tools/compare_irt_geant4.py`).
- `resampling = "incremental"` keeps the sampled reaction times and draws only the pairs of new molecules: 100 overlapping tracks reach 1 us in half a second per trajectory. `"full"` (default) draws every pair again after each reaction.
- Reactions of one molecule in the IRT stage: first order, with water, and with a background species (a species with an initial concentration, such as O2 or the ions of the pH). The background keeps its concentration; what it loses is reported. `background_reactions = false` turns them off.
- Several tracks can be injected at their own times and offsets. `boundary = "reflective"` puts walls at the faces of a box, with mirror images for pairs near a wall.
- With an `[irt]` section, `radkin voxel --solver ssa` runs the IRT stage first and hands its molecules over to the RDME at `irt.t_end`; `yields.csv` covers both stages.
- Reaction entries take an optional `irt_type = "full" | "partial"` (fully or partially diffusion-controlled). Without it the type comes from the rate constant against the diffusion limit. Species radii and charges can be set in `[irt.species.<label>]`.
- New `data/chem_g4.toml`: the Geant4-DNA 11.4.3 option3 chemistry, written by `tools/export_g4_chem.py` from the reaction table Geant4 dumps. Its self-reactions (`2 OH_r -> H2O2` and five others) hold twice Geant4's rate constant, as radkin's convention for `2 A` needs; the export of the earlier comparison sessions ran them at half Geant4's rate.
