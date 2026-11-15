# Session: geant4-single-track-gvalues-10us

**Date:** 2026-11-15T00:24:40.100Z
**Status:** Grill phase complete

## Problem Statement

Test radkin's newest voxel model, the RDME (`solver = "ssa"`, mesh-corrected rates, species levels), on single Geant4-DNA tracks. It should give G(t) up to 10 µs, starting from five handover times, with an initial voxel size guessed and then checked. The deterministic tests before this (sessions 2026-11-13 and 2026-11-14) could not resolve spur recombination: at 40 nm, seeded at 1 ps, G(e_aq) was 3.17 at 10 µs against 2.63 for IRT. The RDME is valid down to the critical spacing (1–4 nm for radical reactions), so it may capture the spur stage. This prepares the next session, which raises the dose rate and checks its impact on G-values.

## Context & Constraints

- **Current behavior:** The Geant4 project `${DEV_DIR}/GEANT4/SIM/radkin` (copy of geant4-rs, 11.4.3) has a 150 µm water cube, or a hard-coded 2×2×10 µm segment in vacuum (`/rs/det/segment`). It uses physics option4, chemistry option3 and IRT_syn. It writes radkin seed CSVs (molecules), edep and G(t) counts per event, and has 10 ns and 50 ns seed macros. In radkin, the RDME stores dense i64 counts per voxel and species (about 180 B per voxel per trajectory) and runs its trajectories in parallel with rayon.
- **Pain point:** Mean-field voxels of 40 nm and up dilute the spurs. The fine voxels the RDME needs only fit in memory (32 GB, 16 cores) on a small grid: 2×2×2 µm at 10 nm is 8e6 voxels, at 5 nm 6.4e7 voxels (about 11 GB, so 2 trajectories at a time).
- **Dependencies:** Geant4 changes stay in the Geant4 project, outside the repo. radkin code is committed only if a fix turns out to be needed; it then follows AGENTS.md (CI clean, golden tolerances unchanged, rsoda pinned by rev). Geant4 option3 IRT cannot scavenge O2 (`G4DNAScavengerMaterial` is used only by the UHDR/Gillespie path), so the aerated cases have no Geant4 reference.
- **Tech stack:** Geant4-DNA C++ (CMake, `run.ps1`); radkin `radkin voxel` with `[ssa]`; Python `-I` scripts to convert seeds, write configs, run sweeps and plot.

## Success Metrics

- Geant4 run on a 2×2×2 µm water segment in vacuum: 100 keV e⁻, IRT_syn, 100 events, the same random seeds in every pass. Five seed passes (handover at 1 ps, 100 ps, 1 ns, 10 ns, 100 ns) plus one G(t) pass to 10 µs (deaerated IRT reference). The edep of event i is identical in every pass.
- A pilot (1 event, 1 ps, 10 nm, deaerated) measures wall time and memory. The ensemble or the sweep is resized if it does not fit.
- RDME G(t) up to 10 µs, events 0–9 × 4 trajectories, G pooled as sum(N)/sum(E/100 eV), net of the background, with SEM:
  - 10 nm: the 5 handovers × 3 chemistries (chem_g4 deaerated, chem_g4 + 52 µM O2, CONFIG_PAPER_01 aerated);
  - 20 nm and 5 nm: 1 ps, deaerated only.
- Plots and `report.md` in `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/`. They compare the RDME with the IRT reference and with the earlier deterministic results (11-13, 11-14), assess voxel convergence and the effect of the handover time, and end with recommendations for the dose-rate session.

## Architecture & Approach

1. **Geant4:** make the segment length configurable, e.g. `/rs/det/segmentZ`, or a segment-size command. The gun starts at the segment's −z face. Write macros for the 5 seed passes and the G(t) pass (common.mac with the same threads and seeds), run them with `run.ps1`, and check that edep matches across passes.
2. **Configs:**
   - `chem_g4.toml` (Geant4-exported, 11-14 session; deaerated, pH 7);
   - `chem_g4` plus 52 µM O2 as an initial concentration;
   - CONFIG_PAPER_01 aerated (52 µM O2, plus the O(³P) reactions).
   - In the aerated runs O2 goes on a coarse species level, taken from radkin's logged suggestion and checked against level 0 on the pilot.
   - Every case has `solver = "ssa"`, the default hellander mesh correction, and a grid equal to the water box.
3. **Runs:** a Python sweep driver over (chemistry, handover, voxel size, event), one `radkin voxel` per event with 4 trajectories, and RAYON_NUM_THREADS set from the memory budget. It launches detached, with file redirection.
4. **Analysis:** pool over events, subtract the background, plot G(t) per species with the IRT reference, and write the summary CSV and the report.

## Assumptions & Trade-offs

- Assumes that reflective (zero-flux) walls at the water-box faces are close enough to the IRT molecules, which diffuse freely. The track runs about 1 µm from the lateral walls, and the e_aq diffusion length at 10 µs is about 0.5 µm.
- Initial voxel guess 10 nm: above h* (1–4 nm), and √(6Dt) for e_aq reaches 10 nm at about 3 ns. 20 and 5 nm check convergence at 1 ps only.
- Deterministic (Strang) runs are not repeated; earlier sessions supply them. No aerated Geant4 reference.
- The 2×2×2 µm segment is a new geometry, so its IRT reference is new. The 11-13 numbers (2×2×10 µm) are compared qualitatively only.

## Open Questions

None. The pilot may resize the ensemble or the 5 nm case.

## Notes

- Earlier results: `.scratch/tests/2026-11-13__geant4-single-track-gvalues/report.md` (deterministic, 1 ps seeds) and `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/report.md` (IRT_syn handover at 10/50 ns).
- RDME docs: `docs/ssa.md`; ADRs 0003–0005.
- Glossary: Handover time (new).
