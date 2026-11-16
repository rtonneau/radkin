# Session: irt-syn-port

**Date:** 2026-11-15T16:32:11.681Z
**Status:** Grill phase complete

## Problem Statement

radkin gets its own IRT solver (new `src/irt/`), starting from the 1 ps molecules, so that Geant4 is only needed up to the end of the physico-chemical stage. Promoted from the Idea `irt-syn-port` (`.work/idea/irt-syn-port.md`).

Why:

- a study needs six Geant4 passes today (five seed passes plus one G(t) pass);
- independent per-track IRT stops being valid above about 1e7 Gy/s, where tracks overlap before the handover;
- Geant4's IRT only samples track-track pairs, so it has no aerated spur stage (no O2 scavenging, no reactions with the background H2O, H3O+, OH-);
- the RDME cannot reproduce spur •OH chemistry from 1 ps by refining the grid (PR #18 report: •OH +40 % at 10 µs, H2O2 0.24 against 0.66).

## Context & Constraints

- **Current behavior:** radkin starts from Geant4-DNA molecules at a handover time chosen in the Geant4 run. Seeds are parsed as molecules with positions (`src/voxel/seed.rs`) and binned onto the voxel grid, where the RDME or the Strang driver computes the chemistry. The 1 ps seed pass matches Geant4 (G(e_aq) 4.422).
- **Pain point:** the particle chemistry between 1 ps and the handover is done by Geant4's IRT_syn, one track at a time, in free space, without background reactions. Evidence: `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/report.md`.
- **Dependencies:**
  - The Geant4 algorithm it starts from (11.4.3: `G4DNAIndependentReactionTimeStepper.cc`, `G4DiffusionControlledReactionModel.cc`, `G4DNAMakeReaction.cc`, `G4DNAMolecularReactionTable.cc`, `G4IRTUtils.cc`):
    1. pair search within a cutoff radius (`G4IRTUtils::GetRCutOff`);
    2. reaction time per pair from the Smoluchowski survival probability: probability `Reff/r`, `t = (r - Reff)^2 / (4 D erfcinv(U/W)^2)`; a partially diffusion-controlled reaction adds a second time built from `kact/(kact + kdif)`;
    3. a pair inside `Reff` at the start time reacts at once with the reaction's probability;
    4. jump to the earliest reaction, diffuse every molecule over that step, place reactants and products at sampled positions;
    5. clear the reaction set and resample every pair after each reaction step.
  - Radii in Geant4: a fully diffusion-controlled reaction has `Reff = kobs / (4 pi D NA)`; a partially diffusion-controlled one uses the sum of the two van der Waals radii, and the Onsager radius from the charges for ion pairs.
  - Geant4's IRT sources have no boundary handling: the wall rule is radkin's own.
  - The Geant4 network and diffusion coefficients (`chem_g4.toml`, `chem_g4.diffusion.json`) and their exporter `export_g4_chem.py` exist only under `.scratch/`, untracked. The export holds the reaction type as a comment only, and no radius or charge.
  - The prerequisite `rdme-reversible-mesh-correction` is done (PRs #19 and #20, ADR 0006).
  - The existing voxel solvers, the SSA random streams (`src/ssa/rng.rs`) and the G(t) outputs are reused, not changed in behaviour.
- **Tech stack:** Rust 2024, `thiserror` (`RadkinError`), rayon, the existing `src/voxel/` and `src/ssa/` modules; Geant4-DNA 11.4.3 IRT_syn as the reference; `tools/compare_geant4.py` for the comparison.

## Success Metrics

- Acceptance: in its single-track, deaerated, free-space configuration with full resampling, G(t) of all species from 1 ps to 10 µs agrees with Geant4 IRT_syn on the same 100 events within statistical error.
- Every reaction of the Geant4 network resolves to a type and a reaction radius equal to Geant4's own values; the agreement of the `k_obs >= k_dif` rule with Geant4's types on that network is reported.
- The sampled reaction times (fully and partially diffusion-controlled, ion pairs, contact pairs) reproduce the analytic survival probabilities in unit tests.
- A first-order decay follows its exponential, and a uniform well-mixed box matches the unit-cell SSA or LSODA result.
- A uniform distribution stays uniform between reflective walls; the survival probability of one pair in a small box agrees with a brute-force Brownian walk with reflection; widely spaced tracks injected at scheduled times reproduce the single-track G(t).
- Incremental resampling agrees with full resampling on the same seeds.
- About 100 overlapping tracks (order 1e5 to 1e6 molecules) reach the handover in minutes; the run time against the number of tracks is measured and recorded.
- One run goes from the IRT stage to the RDME at the configured handover time with molecule counts conserved; a track scheduled after the handover is a config error.
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings` and `cargo test --release` stay clean, and no golden tolerance is loosened.

## Architecture & Approach

A new module `src/irt/` holds the IRT stage: every molecule is a particle with a position, and reactions are drawn as independent reaction times per pair, and per molecule for background reactions. It uses Geant4's IRT sampling formulas but is designed for many overlapping tracks from the start. The Geant4 match is the acceptance test, not the design constraint.

Decisions:

1. **Resampling.** Two selectable rules. Full resampling is the reference mode and carries the Geant4 match. Incremental resampling is accepted when it agrees with the full mode on the same seeds.
2. **Networks.** Any radkin network runs in the IRT stage, not only Geant4 exports.
   - A reaction type given by the network (a Geant4 export) wins.
   - A reaction without a type gets one from its rate constant, with no free parameter: `k_dif = 4 pi NA D R` from the species table (R = sum of the radii, Onsager-corrected for ions). `k_obs >= k_dif` means fully diffusion-controlled; lower means partially diffusion-controlled, with `k_act` derived from `k_obs` and `k_dif`.
   - radkin ships one species table (radius, charge, diffusion coefficient) for the water-radiolysis species, taken from Geant4. A network may override an entry.
   - A species missing from the table makes its reactions fall back to fully diffusion-controlled, with a warning.
3. **Background reactions.** First-order reactions with fixed-concentration species, as exponential reaction times per molecule.
4. **Scheduled tracks.** Tracks are injected at scheduled times (a time and an offset per seed, as the voxel solvers accept).
5. **Walls.** Reflective walls, consistent with the RDME box. They reflect diffusion displacements and product placement. A pair within the cutoff radius of a wall also samples a reaction time against the partner's mirror image. Images are taken across every wall within the cutoff radius, including the edge and corner images when two or three walls are that close.
6. **Pair search.** A spatial grid from the start, sized by the cutoff radius, for about 100 overlapping tracks.
7. **Handover.** One fixed handover time per run, set in the config; the handover to the RDME is in-process. A track scheduled to arrive after it is a config error. A long pulse therefore means a late handover or an IRT-only run.
8. **Config and CLI.** An `[irt]` table (resampling rule, handover time, species-table overrides). `radkin irt` runs the IRT stage alone; `radkin voxel` with `[irt]` runs it up to the handover.
9. **Reuse.** Random streams and G(t) outputs reuse the SSA ensemble's (`src/ssa/rng.rs`, same CSV columns).
10. **Repo data.** The Geant4 network, its diffusion file and `export_g4_chem.py` move from `.scratch/` into `data/` and `tools/`.
11. **Cut at 1 ps.** Going earlier means porting the dissociation channels, `G4DNAWaterDissociationDisplacer`, electron-hole recombination and electron thermalisation, which buys nothing for dose rate.

Milestones, all in this session (one branch, one PR):

1. **Species table and reaction typing.** Every reaction of the Geant4 network resolves to a type and a reaction radius equal to Geant4's own values; the rate-constant rule is compared with Geant4's types on that network.
2. **Pair sampling kernel.** The sampled reaction times reproduce the analytic survival probabilities in unit tests.
3. **Single-track solver, full resampling.** The acceptance test against Geant4 IRT_syn.
4. **Background reactions.** Exponential decay and well-mixed box limit cases.
5. **Walls and scheduled tracks.** Uniform stays uniform, the one-pair Brownian-walk check of the image rule, and widely spaced tracks reproduce the single-track G(t).
6. **Incremental resampling.** Agreement with the full mode on the same seeds; run time against the number of overlapping tracks measured and recorded.
7. **Handover to the RDME.** One run from the IRT stage to the RDME with molecule counts conserved; a track scheduled after the handover is a config error.
8. **Docs.** ADR 0007 and `docs/irt.md` describe the solver, its config and the provisional status of multi-track aerated results.

## Assumptions & Trade-offs

- Parts without a Geant4 IRT reference (background reactions, walls, overlapping tracks) are checked on limit cases against radkin's own solvers only. Multi-track aerated results stay provisional until the Idea `geant4-uhdr-multi-track-reference` delivers a comparison.
- IRT treats pairs as independent. This gets worse as tracks overlap densely, and nothing in scope measures that error.
- The image rule treats the image as one more independent pair, the same approximation IRT already makes; the Brownian-walk check measures it for one pair only.
- The `k_obs >= k_dif` rule may disagree with Geant4's own reaction types on the Geant4 network; milestone 1 reports it, and the network's own type always wins.
- Networks that do not come from Geant4 have no reference at the particle level.
- With a single handover time, long pulses force a late handover or an IRT-only run, which may be slow.
- Incremental resampling departs from the Geant4 algorithm; it is only as good as its agreement with the full mode on the cases tested.
- The scope is large for one session: eight milestones behind one acceptance test, in one PR. Chosen knowingly over splitting after milestone 3 or 5.
- Not doing: anything before 1 ps; a target of about 1000 tracks (order 1e7 molecules) with parallel pair search; several handover times per run.

## Open Questions

None blocking. Two results are measured during the work, not decided now:

- how often the `k_obs >= k_dif` rule agrees with Geant4's types (milestone 1);
- the run time against the number of overlapping tracks (milestone 6).

## Notes

- The glossary already holds the terms settled in the Idea's grill: **Handover time** (extended to the in-process IRT stage), **IRT stage**, **Full resampling**, **Incremental resampling** (`.work/GLOSSARY.md`, uncommitted at the start of this session).
- The design gate was passed on 2026-11-15: the Idea's decisions were confirmed as a whole, and four questions were re-grilled (type rule, wall rule, scale, session scope).
- `G4IRTUtils::GetRCutOff` uses a 1000 ns cutoff time, a 1 % reaction probability, a 1.45 nm maximum reaction radius and twice 9.46e9 nm²/s as the maximum relative diffusion coefficient.
