---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

CLI flags for beam and initial concentrations, removing the wrapper TOML in comparisons

- Phase-1 ruling: radkin CLI has no --beam-type/--beam-param/--initial; compare_radiopyo.py writes a wrapper TOML that includes the config.
- Flags must still go through radiopyo's merge rules (BeamCollection merge, default-beam drop) to stay comparable.

## Changelog

- 2026-11-13 — capture — CLI flags for beam and initial concentrations, removing the wrapper TOML in comparisons (radkin @ main)
