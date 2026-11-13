---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

Close the Phase-1 acceptance gaps: long pulsed spans in run mode and equilibrium net yields

- Run mode with pulsed/macro-pulsed beams is only checked on the first pulse(s) with max_step 1e-6; chunked mode covers full spans.
- Equilibrium net yields are compared with a floor of rtol x gross flux, which could hide an error below 1e-6 of that flux; comparing forward and backward gross yields separately would be stricter.
- Only radiopyo-parity configs were exercised; no randomized/property tests of parser or merge rules.

## Changelog

- 2026-11-13 — capture — Close the Phase-1 acceptance gaps: long pulsed spans in run mode and equilibrium net yields (radkin @ main)
