# Ticket 04: strang-adaptive-driver

**Status:** ✅ Done

## Local Test Result

`cargo test --release voxel::driver` → 3 passed; all voxel tests 16 passed; clippy clean.

## Review Notes

Base 7c7fe68. `VoxelRun`: Strang R(dt/2) D(dt) R(dt/2). Step doubling compares one dt step with two dt/2 steps, keeps the two-half-step state, uses the scaled max-norm error, and sets factor = safety·err^(-1/3) clamped to [shrink, grow]. `advance_to` lands exactly on the target, so events are handled by the caller between calls. `inject` adds amounts, sets the Rannacher flag (BE substeps for the next accepted step's diffusion) and resets dt. Tests: a uniform field through the full driver matches one Phase-1 solve at rtol 1e-6 (12 voxels, two advance targets); on an A+B→C 1D interface, the error falls when rtol goes 1e-2 → 1e-4 against a 1e-7 reference; an injection is conserved, spreads, and stays non-negative. Accepted inline.

## Blockers / Challenges

None

## Commits

- 5b1db69 feat(voxel): Strang splitting driver with step-doubling adaptive dt and injections (ticket 04)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
