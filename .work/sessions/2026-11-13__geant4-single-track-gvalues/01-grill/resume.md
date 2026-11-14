# Session: geant4-single-track-gvalues

**Date:** 2026-11-13T22:57:00.828Z
**Status:** Grill phase complete

## Problem Statement

Test radkin's voxel reaction-diffusion (Phase 2) on real Geant4-DNA single tracks: G-value time evolution from 1 ps to 10 µs for 100 keV electrons in the 0.5 mm water box of geant4-rs (`${DEV_DIR}/GEANT4/SIM/geant4-rs`), with an initial voxel size guessed and then swept, against a Geant4-DNA IRT reference. This prepares the follow-up session (higher dose rate: several tracks in the box, impact on G-values).

## Context & Constraints

- **Current behavior:** geant4-rs (not a git repo) shoots 100 keV e⁻ along +z from the -z face of a 0.5 mm water cube (DNA physics option4, chemistry option3, IRT), stops the chemistry at 1.001 ps and writes the molecules (radkin seed CSV, nm) per thread; 100 events, ~12k molecules per event, track extent ~36×37×71 µm (max ~135 µm). No deposited energy is written. radkin `radkin voxel` records N(t) and G(t) at every splitting step in `yields.csv` (G needs `energy_ev` per seed); rsoda runs on every voxel each half-step, empty voxels included.
- **Pain point:** spur recombination happens at nm–tens of nm in 1–100 ns, while diffusion over 10 µs is only ~0.6 µm; a whole track cannot be resolved at spur scale on a uniform grid, and at µm voxels a track voxel holds ~10 nM so G(t) stays nearly flat apart from scavenging. There is no reference G(t) to 10 µs and no edep for G-values.
- **Dependencies:** rsoda pinned by rev (no path dependency committed); radiopyo stays the Phase-1 reference; golden tolerances unchanged; kernel does not allocate; CI (fmt, clippy -D warnings, tests) clean on the session branch with a PR. Geant4 11.4 install from `G4_ROOT`, sources at `${DEV_DIR}/GEANT4/geant4-v<X>-source` for option3 reactions and molecule-counter patterns (chem6, molcounters examples).
- **Tech stack:** Rust 2024 (radkin, rayon), Geant4-DNA C++ (geant4-rs, CMake, `run.ps1`), Python `-I` for conversion and plotting scripts.

## Success Metrics

- geant4-rs writes the deposited energy per event and, with the chemistry run to 10 µs, the Geant4 IRT G(t) on log-spaced times for both setups (whole track; track segment run).
- radkin advances background voxels by copying one shared background integration; active voxels give unchanged results; fmt, clippy and tests clean; a test covers the background-voxel path.
- radkin G(t) from 1 ps to 10 µs for the whole track (event 0) at 2 µm / 1 µm / 500 nm and for the track segment at 80 / 40 / 20 nm, in both chemistries (option3-matched deaerated, CONFIG_PAPER_01 aerated), compared against Geant4 IRT in plots and a short report in `.scratch`.

## Architecture & Approach

1. **Geant4 (geant4-rs):** add edep scoring per event (and the region's edep for the segment setup). Run the IRT chemistry to 10 µs and record G(t) with the Geant4 molecule counter on log-spaced times. Keep writing the 1.001 ps seed molecules (seed and reference must come from the same events, e.g. a seed-only pass and a full-chemistry pass with the same random seeds, or both in one run). Two setups: (a) the existing whole-track run (0.5 mm box); (b) a **track segment run**: water only in a 2×2×10 µm region at the track entrance, vacuum world, ~5 events, so edep, molecules and IRT G(t) all refer to that region.
2. **radkin code:** **background voxels**: a voxel whose state matches the shared background vector within a tolerance far below the solver atol is advanced by copying one shared background integration per half-step instead of its own rsoda run; active voxels are integrated exactly as before.
3. **Chemistry configs:** an option3-matched config (reactions and rates read from `G4EmDNAChemistry_option3.cc`, deaerated, no background species), and CONFIG_PAPER_01 aerated (52 µM O2, with the O(3P) reactions as in `voxel_example.toml`) as the physics case.
4. **Runs:** whole track: event 0, grid = bounding box + 2 µm margin, voxels 2 µm / 1 µm (initial guess) / 500 nm. Segment: grid = region + 1 µm margin, voxels 80 / 40 (initial guess) / 20 nm, ~5 events run separately and averaged. A short benchmark first to predict runtimes; overnight runs are allowed.
5. **Output:** seed/edep conversion, run and plotting scripts, configs, plots and a short report, all under `.scratch/tests/2026-11-13__geant4-single-track-gvalues/`; only the radkin code change is committed.

## Assumptions & Trade-offs

- Mean-field per voxel: µm voxels dilute the spurs (little recombination); only the segment runs approach spur scale. The whole-track runs show the voxel-size dependence, not a converged answer.
- IRT has no walls, radkin has zero-flux walls: the margins (2 µm whole track, 1 µm segment) exceed the ~0.6 µm diffusion length at 10 µs, so the walls barely matter.
- Rate constants: radkin uses bulk rates, Geant4 IRT uses its reaction-type models; small differences are expected even for matched chemistry.
- Not doing: dose rate or multiple tracks (next session), the full 0.5 mm box as grid, non-uniform or adaptive grids, committing study assets (configs, scripts, report stay in `.scratch`).

## Open Questions

- The exact Geant4 molecule-counter API in 11.4 and whether IRT to 10 µs on a full 100 keV track runs in reasonable time: checked against the Geant4 sources and examples during implementation.
- The background-voxel tolerance value and how it interacts with ADI round-off (a uniform field may not stay bit-exact after diffusion): settled in the radkin ticket.
- The actual per-voxel cost, hence which sweep points are practical: settled by the benchmark.

## Notes

- Initial voxel-size guesses (the user asked for them): 1 µm for the whole track, 40 nm for the segment.
- The README's mean-field caveat says large voxels overestimate recombination; dilution argues they underestimate it. Check against the results and correct the README if they confirm it (only if the change is clear).
- Decisions by the user: both strategies compared; Geant4 IRT reference to 10 µs; both chemistries; edep added to Geant4; dedicated segment run; overnight runs fine; background-voxel code change regardless of the benchmark; study assets in `.scratch`; run with auto.
