# Ticket 02: adi-diffusion

**Status:** ✅ Done

## Local Test Result

`cargo test --release voxel::diffusion` → 4 passed (gaussian_3d 61³, 1D/2D degenerate, mass+positivity, second order in time ratio≈4); clippy clean.

## Review Notes

Base 83758df (ticket-start was recorded after implementation because the first start call failed in a broken shell command; the work started from 83758df). ADI is the product of 1D theta-scheme factors per axis (LOD). Ax/Ay/Az commute on a uniform Neumann grid, so CN is 2nd order (verified: error ratio ~4 on halving dt), unconditionally stable and mass-conserving (1e-12). Backward-Euler factors are M-matrix inverses, so positivity holds with the Rannacher start. Gaussian tolerance 1% of peak with σ0=4 nm (at 2 nm the 2.2% 3D error is the O(dx²/σ²) spatial error). No allocation per step. Accepted inline.

## Blockers / Challenges

The first 3D Gaussian test was limited by spatial resolution, fixed with a wider Gaussian.

## Commits

- 212f4c1 feat(voxel): Crank-Nicolson ADI diffusion with zero-flux boundaries and Rannacher start (ticket 02)

## Time Spent

0m (ticket-start.js to ticket-complete.js)
