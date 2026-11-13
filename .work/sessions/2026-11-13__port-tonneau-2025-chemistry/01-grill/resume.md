# Session: port-tonneau-2025-chemistry

**Date:** 2026-11-13T22:41:04.748Z
**Status:** Grill phase complete

## Problem Statement

radkin has no ready-to-run config for the Tonneau 2025 (PMB 70 235021) chemistry. The original-runs chemistry is bundled verbatim as `data/CONFIG_PAPER_01.toml`, but the inputs that paper_runs.py sets in code are not in any config: variable G-values, LET 0.56 keV/µm and a single-pulse beam. HO3, a species of this chemistry, also has no default diffusion coefficient for voxel runs.

## Context & Constraints

- **Current behavior:** `data/` matches radiopyo's TOMLs byte for byte (line endings aside). Every run needs a config path. `DEFAULT_DIFFUSION` (src/voxel/config.rs) has 15 Geant4-DNA option3 species but no HO3. Geant4-DNA 11.4.3 has no HO3 at all.
- **Pain point:** Reproducing the paper's original runs needs the user to set variable_GValues, the LET and the beam by hand. HO3 does not diffuse in voxel runs.
- **Dependencies:** Keep CONFIG_PAPER_01.toml unchanged as the single definition of the chemistry. Keep the golden tolerances.
- **Tech stack:** Rust, TOML config with `[includes]` (the `radkin/` prefix points to the bundled data folder).

## Success Metrics

- `data/tonneau_2025.toml` loads with `UnitCell::from_toml` and `radkin` runs it without extra flags.
- Its settings match paper_runs scenario "original": the CONFIG_PAPER_01 chemistry, variable_GValues = true, LET 0.56, and a single-pulse beam of 30 Gy at 1e7 Gy/s.
- HO3 gets 2.0e-9 m²/s by default in `[voxel]` runs. A test checks this.
- fmt, clippy and the tests are clean.

## Architecture & Approach

A thin `data/tonneau_2025.toml` with `[includes] toml = ["radkin/CONFIG_PAPER_01.toml"]`, `[parameters] variable_GValues = true` and one `[[beam.single_pulse]]` (dose rate 1e7 Gy/s, max dose 30 Gy, LET 0.56). Add `("HO3", 2.0e-9)` to `DEFAULT_DIFFUSION` with a comment: Geant4-DNA has no HO3, so it takes O3's value. Update README and data/README.

## Assumptions & Trade-offs

The variant is the "original runs" one (CONFIG_PAPER_01, LET 0.56), not the paper-methods one (pure-water config, LET 1.17). There is no embedded default config and no port of the paper grid. HO3's coefficient is an assumption, not a measured value. The default beam is one point of the paper grid; others are chosen with CLI flags or an own TOML.

## Open Questions

None.

## Notes

The user chose: a bundled TOML only, the original-runs variant, HO3 = 2.0e-9 like O3, and a thin include file.
