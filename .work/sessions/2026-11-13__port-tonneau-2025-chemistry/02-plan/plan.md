# Implementation Plan

**Session:** port-tonneau-2025-chemistry
**Date:** 2026-11-13T22:46:44.219Z
**Estimated effort:** 1 hour

## Strategy

First add the bundled Tonneau 2025 config and its test, then the HO3 diffusion default. The two changes are independent and small.

## Tickets Overview

- **Ticket 1:** Add data/tonneau_2025.toml (includes CONFIG_PAPER_01, with variable G-values, LET 0.56 and a single pulse of 30 Gy at 1e7 Gy/s), plus a test and docs.
- **Ticket 2:** Add HO3 at 2.0e-9 m²/s to DEFAULT_DIFFUSION, with a test and docs.

## Sequencing Rationale

They don't depend on each other. The config comes first because it is the main deliverable.

## Risks & Mitigation

- **Risk:** the beam TOML key names differ from radiopyo's → **Mitigation:** copy the keys of an existing single-pulse config or test.

## Assumptions

- The `radkin/` include prefix resolves to the bundled data folder.
