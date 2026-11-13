# Ticket 03: voxel-reaction-step

**Status:** ✅ Done

## Local Test Result

`cargo test --release voxel::reaction` and `--no-default-features` → 2 passed each; full suite ok; clippy clean with and without the `parallel` feature.

## Review Notes

Base baa1c65. `reaction_step` gathers voxels to AoS, solves each with `solver::solve` (no dose source, t_eval = end) and scatters back. rayon `par_chunks_mut` behind the default `parallel` feature. The test asserts bit equality with one-voxel-at-a-time solves. A uniform 8-voxel grid stepped over 12 log-spaced reaction steps matches one Phase-1 solve over (0, 1e-3) at rtol 1e-6 (CONFIG_PAPER_01 chemistry, e_aq decays >50%, asserted). VoxelState units are µmol/L like the kernel. No rsoda change needed: a fresh Lsoda per voxel per step. That restart cost is the expected hotspot and is left for a benchmark.

## Blockers / Challenges

Windows Python wrote µ as cp1252; fixed. Closure needed Send for rayon.

## Commits

- 9847303 feat(voxel): per-voxel rsoda reaction step, rayon-parallel behind default feature (ticket 03)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
