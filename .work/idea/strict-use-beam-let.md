---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

Make use_beam(label, None) match radiopyo (error) or document the LET fallback

- Phase-1 ruling: use_beam without a LET falls back to the beam's own LET, where radiopyo raises.
- Cost: a run proceeds where radiopyo would stop; decide whether parity or convenience wins.

## Changelog

- 2026-11-13 — capture — Make use_beam(label, None) match radiopyo (error) or document the LET fallback (radkin @ main)
