# Implementation Plan

**Session:** rdme-mesh-corrected-rates
**Date:** 2026-11-14T23:27:19.360Z
**Estimated effort:** 1 day

## Strategy

Add the config key first, then the correction itself in `RdmeSetup` (per-level factor, pair constant shared with `validity_bound`, fallback and warning rules in `run_rdme`), then validate it against the association-time model problem and update the docs with the rerun small-voxel table.

## Tickets Overview

- **Ticket 1:** `[ssa] rate_correction = "hellander" | "none"` parsed (unset = default), and an explicit value refused on unit-cell runs.
- **Ticket 2:** Mesh correction applied per coarsening level in `RdmeSetup`, with `kappa = 2k` for `2 A` (also in `validity_bound`), `h <= h*` left uncorrected with a warning, and the anisotropic fallback/error.
- **Ticket 3:** One-A-one-B association-time validation test, and `docs/ssa.md` (key, critical spacing, small-voxel table rerun with and without the correction).

## Sequencing Rationale

Ticket 2 needs the parsed setting from ticket 1. Ticket 3 validates and documents the behaviour of ticket 2 and reruns the example on it.

## Risks & Mitigation

- **Risk:** Default on changes existing RDME integration tests' results → **Mitigation:** their voxels are coarse (factor near 1); if one fails, check whether the change is the expected physics before touching the test, never loosen golden tolerances.
- **Risk:** The association-time test is slow or noisy → **Mitigation:** small grid (e.g. 8^3 to 12^3), `D_A = 0`, many cheap trajectories; tolerance set from the measured SEM, and the uncorrected run must sit several SEM away.
- **Risk:** Zero-flux boundaries differ from the periodic domain of Theorem 1 → **Mitigation:** pick `L >> h`, compare with `V/k` plus the `tau_D` term where needed, and record the observed agreement in the test comment.

## Assumptions

- Each `Level` keeps its own `SsaModel`, so the correction can be baked into `c` at setup; tau-leaping and coarsening read `c` unchanged.
- `RdmeSetup::new` is crate-internal API used only by `run_rdme` and tests, so its signature can grow a correction argument.
