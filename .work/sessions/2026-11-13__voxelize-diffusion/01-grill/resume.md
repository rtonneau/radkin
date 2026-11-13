# Session: voxelize-diffusion

**Date:** 2026-11-13T21:01:03.926Z
**Status:** Grill phase complete

## Problem Statement

radkin only simulates homogeneous chemistry (one unit cell). The goal is to simulate the physical track of electrons or protons in Geant4-DNA, then use its output to seed a voxelised simulation whose voxels are small enough to resolve heterogeneous track chemistry: solve ODE, diffuse, solve ODE, and so on.

## Context & Constraints

- **Current behavior:** `UnitCell` integrates one well-mixed state vector with rsoda. Species are created by a Beam through G-values. The kernel (`src/kernel.rs`, `src/matrices.rs`) already works on one cell at a time with buffers passed in by the caller.
- **Pain point:** there is no spatial dimension, no diffusion, and no way to import Geant4-DNA species distributions.
- **Dependencies:** keep the Phase-1 results and golden tests unchanged. Reuse the per-cell kernel unchanged. rsoda stays a pinned git dependency. ADR 0001 (pluggable linear solver) stays available for a later monolithic reference solver.
- **Tech stack:** Rust 2024, rsoda (LSODA), thiserror, TOML config, rayon (new, behind a default-on cargo feature), CSV.

## Success Metrics

- A spatially uniform seed on the voxel grid reproduces the Phase-1 unit-cell result within rtol 1e-6 (per voxel and integrated).
- Diffusion alone (no reactions) of a point source matches the analytic 3D Gaussian within a stated tolerance, and total mass is conserved to round-off with zero-flux boundaries.
- Parallel (rayon) and serial runs are bit-identical.
- An end-to-end run from a Geant4-DNA molecule CSV produces snapshots and domain-integrated G(t) CSV. Its G(t) is compared qualitatively with Geant4-DNA chem in a documented tool/notebook.
- `cargo fmt --check`, `cargo clippy --release --all-targets -D warnings` and `cargo test --release` stay clean on Ubuntu and Windows.

## Architecture & Approach

- **Voxel grid:** uniform Cartesian nx×ny×nz with per-axis spacing and origin, configured in TOML. Each voxel holds one full species state vector (structure-of-arrays storage).
- **Seed import:** a CSV of molecules (species, x, y, z in nm, optional event_id). Names are mapped to radkin species through a TOML `[geant4.aliases]` table, with a default table for the standard Geant4-DNA water species. Unknown names are an error unless listed as ignored. Each molecule is binned into its voxel as 1/(N_A·V). Molecules outside the grid are an error that reports the count and bounding box, unless `clip = true`, which drops them and reports the lost fraction.
- **Seed schedule:** a `[[seed]]` list in TOML, each entry giving a file, an optional event_id filter, an injection time and an optional spatial offset. Homogeneous background species (O2, pH, ...) come from the existing config.
- **Diffusion coefficients:** an optional per-species field (m²/s) in the TOML. Default values match Geant4-DNA. Species without a value do not diffuse.
- **Time stepping:** Strang splitting, reaction(dt/2) → diffusion(dt) → reaction(dt/2). The reaction half-step runs rsoda per voxel, in parallel through rayon, with bit-identical results. Diffusion is Crank-Nicolson ADI (Douglas) with zero-flux Neumann boundaries, one tridiagonal (Thomas) solve per line. A few backward-Euler steps after each injection damp CN oscillations on spiky data (Rannacher start). dt is adaptive by step-doubling (one full step against two half steps), with the step-size controller clamped to land exactly on injection and snapshot times.
- **Output:** snapshots (3D concentration fields per species at requested times, in a simple documented binary format such as .npy) plus a CSV of domain-integrated amounts and G(t).
- **Validation:** uniform seed against the Phase-1 run, diffusion against the analytic Gaussian with mass conservation, and a G(t) comparison against Geant4-DNA chem.

## Assumptions & Trade-offs

- The model is mean-field and continuum. At nm voxels one molecule is ~mM, and a lone molecule can react with itself through second-order terms. The choice of voxel size must be documented and checked; this is a known limitation, not something to fix in Phase 2.
- Splitting error is accepted in exchange for linear scaling and reuse of the kernel. A monolithic method-of-lines solve is not built now (ADR 0002).
- Only a uniform Cartesian grid with zero-flux boundaries: no cylindrical, adaptive, Dirichlet or periodic grids.
- No ROOT reader. The Geant4 side writes the CSV.
- No homogeneous beam source during voxel runs. Sources are seeds only.

## Open Questions

- Which molecule time the Geant4-DNA exporter uses (end of prechemistry ~1 ps vs a later chem time). The seed's injection time carries it, so radkin is agnostic.
- Snapshot file format details (.npy per species/time vs one file per time) are an implementation choice for the plan.

## Notes

- ADR 0002 records the splitting decision. Glossary terms added: Voxel grid, Seed, Splitting step, Snapshot.
- Supersedes the idea `.work/idea/phase2-voxel-diffusion.md` (which leaned to "one ODE system").
