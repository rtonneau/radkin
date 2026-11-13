---
status: captured
scope: project
created: 2026-11-13
reviews: 0
promoted_to: null
---

## Problem

rsoda: sparse or Krylov linear algebra so large stiff systems (Phase-2 voxels) are tractable

- rsoda only has dense and banded Jacobians; n_voxels x n_species unknowns make dense/banded LU unusable.
- Alternative is operator splitting (rsoda per voxel + diffusion step); this Idea is the monolithic-solver path.
- Candidates: sparse LU (e.g. faer/KLU-style) via a user-supplied linear-solver trait, or GMRES with preconditioner.

## Changelog

- 2026-11-13 — capture — rsoda: sparse or Krylov linear algebra so large stiff systems (Phase-2 voxels) are tractable (radkin @ main)
