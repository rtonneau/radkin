# Ticket 01: background-voxels

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`: no fmt diff, no clippy warning, every test binary ok (lib: 47 passed, including `background_voxels_match_full_integration`, `unseeded_grid_is_all_background`, `diffusion_keeps_a_uniform_background_exact`, `background_voxels_setting`).

Extra check: `radkin voxel data/voxel_example.toml` with and without `background_voxels`: same 48 yield rows, max relative difference 9.5e-10.

## Review Notes

Base: 69ee0d8b8477059a9617684c54d76df34e549aa6 (clean tree apart from the gps session record). Inline implementation, no delegate: a localized change in `src/voxel/`, assessed directly from the diff and the tests.

- Shared background: `reaction_step_background` integrates the background vector once per half-step, always, so it stays current. Voxels that match it within `1e-3*atol + 1e-12*|b|` get the result instead of their own rsoda run. The step-doubling branches (big, small) each carry their own copy, and the accepted step keeps the small one.
- Active voxels go through the same `solve` call with the same options. Without a background, the code path is the old one (`reaction_step` is a wrapper), and the existing voxel tests pass unchanged.
- Deviation from the plan: ADI round-off moved a uniform field by about 1e-11 relative over 40 steps, more than the 1e-12 tolerance. With a background set, diffusion is now applied to `y - bg`, which is the same thing mathematically for zero-flux diffusion, and background voxels stay exact (the test asserts equality). So active voxels are not bit-identical to the old code when the background is on: they differ by round-off from the shift, and the solver part is unchanged.
- Counters `background_voxel_steps` and `voxel_steps` are on `VoxelRun`, and the background fraction is logged at info level. Config key `[voxel.run] background_voxels` (default true). README updated.

## Blockers / Challenges

Uniform-field round-off in the ADI solve (see Review Notes), solved by diffusing the deviation from the background.

## Commits

- 1d22775 perf: Share one integration for background voxels (ticket 01)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
