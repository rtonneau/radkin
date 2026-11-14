# Session: geant4-radkin-irt-syn-gvalues

**Date:** 2026-11-14T07:23:44.655Z
**Status:** Grill phase complete

## Problem Statement

Start radkin's voxel chemistry from a later Geant4-DNA handover. A new Geant4 project dedicated to radkin (`${DEV_DIR}/GEANT4/SIM/radkin`) runs IRT_syn to 10 ns and 50 ns and writes the molecules as radkin seeds. Geant4's chemistry is exported as a radkin TOML. radkin then computes G(t) up to 10 µs for single 100 keV electron tracks, with an initial voxel size guessed and then swept. This prepares the next session: higher dose rate and its impact on G-values.

## Context & Constraints

- **Current behavior:** geant4-rs (0.5 mm water cube, 100 keV e⁻ along +z, DNA physics option4, chemistry option3, plain IRT) writes seeds at 1.001 ps plus edep per event. The previous session (`2026-11-13__geant4-single-track-gvalues`) compared radkin started at 1 ps against IRT G(t) to 10 µs. Seeded at 1 ps, radkin stays well above IRT at every voxel size, because spur recombination (nm, 1–100 ns) is not resolved. Its 20 nm sweep was stopped by the user, and that session is left unfinished.
- **Pain point:** a 1 ps handover puts unresolved spur recombination into mean-field voxels. Handing over after most spur chemistry (10–50 ns, spurs about 25–40 nm) should let µm-scale voxels work. The hand-written option3 TOML can also drift from what Geant4 actually uses.
- **Dependencies:** rsoda is pinned by rev, with no path dependency committed. Golden tolerances are unchanged, the kernel does not allocate, and CI (fmt, clippy -D warnings, tests) stays clean. Geant4 11.4 is installed at `G4_ROOT`, with sources at `${DEV_DIR}/GEANT4/geant4-v<X>-source` (IRT_syn, reaction table API). The work goes directly on the current branch `perf/background-voxels` (user).
- **Tech stack:** Rust 2024 (radkin, rayon), Geant4-DNA C++ (CMake, `run.ps1`), Python `-I` for conversion, export and plotting scripts.

## Success Metrics

- `${DEV_DIR}/GEANT4/SIM/radkin` builds and runs: a copy of geant4-rs with a 150 µm water cube, `TimeStepModel IRT_syn`, and 100 events. It writes radkin seed CSVs plus edep per event at 10 ns and at 50 ns, from the same events.
- A Geant4 reaction-table dump (species, diffusion coefficients, reactions, rates, reaction types) is converted by a script into a radkin TOML that radkin loads without error.
- radkin G(t) up to 10 µs for events 0–4 (pooled), at 250 nm (initial guess), then 500 nm and 125 nm. Each set runs for both handover times and both chemistries: the Geant4-exported config (deaerated) and CONFIG_PAPER_01 aerated (52 µM O₂ plus O(³P)). Plots and a short report in `.scratch` compare the results with the previous session's IRT G(t).

## Architecture & Approach

1. **Geant4 project:** copy geant4-rs to `${DEV_DIR}/GEANT4/SIM/radkin`, then make these changes:
   - water cube shrunk to 150 µm;
   - TimeStepModel set to IRT_syn;
   - chemistry stopped and molecules written at 10 ns and 50 ns: two passes with the same random seeds, or one run with two snapshots if the API allows it;
   - edep per event kept;
   - physics and the other settings match geant4-rs.
2. **Chemistry export:** at chemistry initialisation, dump G4DNAMolecularReactionTable and the molecule definitions (D, charge) to a file. A script maps them to radkin TOML species and reactions (rates in radkin units, deaerated).
3. **radkin runs:**
   - events 0–4, each on its bounding box plus a 2 µm margin (not the full 150 µm box: 216M voxels at 250 nm);
   - from the handover time to 10 µs;
   - voxel sweep 250 → 500 / 125 nm;
   - G pooled as sum(N) / sum(E/100 eV).
4. **Output:** scripts, configs, plots and report under `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/`. The Geant4 project lives outside the repo. Only radkin code fixes, if any are needed, are committed.

## Assumptions & Trade-offs

- The reference is the old plain-IRT G(t) (0.5 mm box, other events). The comparison is pooled G per 100 eV, not event by event. IRT and IRT_syn differ slightly.
- In a 150 µm box some tracks escape. Edep scoring keeps G correct.
- Mean field per voxel: 250 nm voxels still dilute the 50 ns spurs a little. The sweep shows how much.
- Not doing: the dose rate (next session), a new Geant4 G(t) reference to 10 µs, the segment geometry, finishing the previous session's 20 nm sweep, committing study assets.

## Open Questions

- Whether IRT_syn can stop at 10/50 ns and give the molecule positions (G4 11.4 API, scheduler end time, ITTrackHolder access): checked in the Geant4 sources.
- How to map Geant4 reaction types (partially/fully diffusion-controlled, with scavengers or bulk (B) species) onto radkin bulk rates: settled in the export ticket.
- The cost of 125 nm (about 50M voxels for a whole track): run only if the benchmark allows it.

## Notes

- Initial voxel-size guess: 250 nm. By 10–50 ns spurs have grown to about 25–40 nm, and 250 nm cost about 12 min per run on the old setup.
- Decisions by the user: copy geant4-rs; seeds at both 10 and 50 ns; TOML generated from Geant4; previous session killed (runs stopped, files kept); 150 µm box; old IRT G(t) reused; 250 nm first; both chemistries; 100 events, radkin on events 0–4; current branch; run with auto.
