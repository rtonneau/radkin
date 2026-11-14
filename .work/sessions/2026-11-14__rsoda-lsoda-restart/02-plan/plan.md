# Implementation Plan

**Session:** rsoda-lsoda-restart
**Date:** 2026-11-14T17:20:06.888Z
**Estimated effort:** 1-2 days

## Strategy

Land the rsoda API first (full-reset `restart`, then the opt-in warm start), each pushed on an rsoda branch so radkin can pin a revision. Then radkin bumps the `rev` and reuses one `Lsoda` per rayon worker in the voxel reaction step, measured under the keep rule. Last, wire the warm start as a `[voxel]` option, off by default, and measure it against the reference run.

## Tickets Overview

- **Ticket 1:** rsoda `Lsoda::restart(t0, y0, h0)` with a full reset equivalent to `new()`, tested in the dense and sparse modes.
- **Ticket 2:** rsoda opt-in warm start on `restart` (keeps the method and the order, falling back to the method only if order > 1 is not valid).
- **Ticket 3:** radkin bumps the rsoda `rev` and reuses one `ChemSystem` and `Lsoda` per worker in `reaction_step_background`; voxel bench under the keep rule.
- **Ticket 4:** radkin `[voxel]` warm-start option (off by default); measure nfev and wall time at equal accuracy, keep or drop.

## Sequencing Rationale

radkin can only pin an rsoda revision that exists on GitHub, so the rsoda tickets come first. Ticket 3 measures the plain restart on its own before the warm start (ticket 4) changes the numerics, so each gain is attributed separately.

## Risks & Mitigation

- **Risk:** an incomplete reset leaks state between voxels → **Mitigation:** a test comparing `restart` with `new()` on several states, in both Jacobian modes, bitwise.
- **Risk:** the gain is within bench noise → **Mitigation:** repeated runs at 1 and 16 threads under the `docs/perf/` keep rule; revert the radkin change if it fails.
- **Risk:** warm start at order > 1 is invalid without history → **Mitigation:** fall back to keeping only the method; drop ticket 4's option if it doesn't win.

## Assumptions

- The rsoda checkout at `${DEV_DIR}/RUST/rsoda` is currently on `perf/per-step-overhead`. The rsoda work goes on a new branch `feat/lsoda-restart`, based on the revision radkin pins (`f5d1312`) or a descendant of it. That branch is pushed to origin so radkin can pin it; nothing is pushed to rsoda `main`.
- A temporary `[patch]` to the local rsoda is used only while developing and is never committed.
