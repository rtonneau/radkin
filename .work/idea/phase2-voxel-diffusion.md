---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

Phase 2: voxelised 3D geometry with species diffusion coefficients solved as one ODE system

- Original project goal: DSMC-like spatial diffusion without individual trajectories, integrated with rsoda.
- Kernel is already per-cell with caller buffers; open choice is operator splitting vs monolithic sparse solve (see [[rsoda-sparse-jacobian]]).
- Needs: geometry/voxel grid, diffusion coefficients per species in TOML, boundary conditions, dose deposition map, outputs per voxel.

## Changelog

- 2026-11-13 — capture — Phase 2: voxelised 3D geometry with species diffusion coefficients solved as one ODE system (radkin @ main)
