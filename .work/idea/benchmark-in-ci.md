---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

Commit the radkin vs radiopyo benchmark and track radkin runtime regressions

- examples/bench.rs and tools/benchmark.py exist but are uncommitted on main.
- radiopyo is not available in CI; a radkin-only timing (criterion or the bench example) could guard regressions.

## Changelog

- 2026-11-13 — capture — Commit the radkin vs radiopyo benchmark and track radkin runtime regressions (radkin @ main)
