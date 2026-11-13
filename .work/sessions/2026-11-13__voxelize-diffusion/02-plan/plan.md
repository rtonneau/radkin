# Implementation Plan

**Session:** voxelize-diffusion
**Date:** 2026-11-13T21:16:35.035Z
**Estimated effort:** 5-8 days

## Strategy

Build bottom-up: first the pieces that can be checked against analytic results (the grid, then pure diffusion), then the per-voxel reaction step, which must match Phase 1. Strang coupling with adaptive dt comes next, followed by the Geant4 seed import and schedule, and finally the outputs, CLI and comparison tooling. Phase-1 code paths and golden tests stay untouched throughout.

## Tickets Overview

- **Ticket 1:** Voxel grid type, grid TOML config and per-species diffusion coefficients.
- **Ticket 2:** Crank-Nicolson ADI diffusion with zero-flux boundaries and a Rannacher start, validated against the analytic Gaussian and mass conservation.
- **Ticket 3:** Per-voxel reaction half-step through rsoda, serial and rayon (bit-identical); a uniform grid matches Phase 1.
- **Ticket 4:** Strang splitting driver with step-doubling adaptive dt aligned to event times.
- **Ticket 5:** Geant4-DNA molecule CSV import, alias table, binning, out-of-grid policy and the `[[seed]]` schedule.
- **Ticket 6:** Snapshot and G(t) outputs, `radkin voxel` CLI subcommand, example config, docs and a Geant4-DNA comparison tool.

## Sequencing Rationale

The grid (1) is needed by everything else. Diffusion (2) and reactions (3) are independent operators, each validated on its own before coupling (4). Seeds (5) need the driver so they can be injected mid-run, and outputs (6) need a full run.

## Risks & Mitigation

- **Risk:** CN oscillations / negative concentrations from single-molecule spikes → **Mitigation:** Rannacher backward-Euler start after each injection; a test asserting non-negativity on a delta seed.
- **Risk:** rsoda per-voxel cost is dominated by startup at tiny dt → **Mitigation:** reuse solver workspaces per voxel and benchmark in `examples/`; skip voxels with a negligible rate (only if exact).
- **Risk:** step-doubling triples the cost → **Mitigation:** reuse the half-step results as the accepted solution (Richardson-free); make tolerances configurable.
- **Risk:** mean-field artefacts at nm voxels → **Mitigation:** document the limitation and voxel-size guidance in the README.

## Assumptions

- rsoda needs no API change; a reset or re-init per voxel is possible through the current API (checked in ticket 3; if not, an rsoda change and a rev bump).
- rayon is acceptable as a new default-on feature dependency.
