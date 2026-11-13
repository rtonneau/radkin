---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

rsoda: expose a public dense-output (DINTDY/dky) API so radkin's t_eval interpolation does not go through integrate(te)

- Phase-1 ruling: t_eval points are interpolated by calling Lsoda::integrate(te) inside the last step instead of a dedicated dky call.
- Cost today: rounding-level differences vs scipy's interpolation and a less explicit contract (relies on integrate not stepping).
- Work lands in rsoda (rtonneau/rsoda), then radkin bumps the pinned rev ae76eee.

## Changelog

- 2026-11-13 — capture — rsoda: expose a public dense-output (DINTDY/dky) API so radkin's t_eval interpolation does not go through integrate(te) (radkin @ main)
