# Ticket 03: association-validation-and-docs

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check` clean; `cargo clippy --release --all-targets -- -D warnings` clean; `cargo test --release`: every suite passes (lib 126 incl. `ssa::rdme::tests::mesh_correction_matches_the_association_time`, 4.0 s).

Small-voxel rerun (`data/voxel_example.toml` copies with `[ssa] trajectories = 64, seed = 1`, outputs in `.scratch/tests/2026-11-15__rdme-mesh-corrected-rates/03-out_*`): `"none"` reproduces the previous table exactly (1.13, 2.53, 2.57, 1.77, 0.68). The corrected run gives e_aq 0.92, OH_r 2.18, H3O_plus 2.06, H2O2 1.75, H2 0.64, against a deterministic 0.84, 2.10, 2.02, 2.01, 0.80.

## Review Notes

Base: ff3e9288a0ba0a60e69338e99478deadb175e35e.

- Deviation from Files to Touch: the association test is a unit test in `src/ssa/rdme.rs`, not in `tests/ssa_rdme.rs`, because it needs the crate-private single-event `Rdme::step` to stop at the first reaction. `cargo test --release ssa` still runs it.
- Probe (4000 trajectories, 11^3 grid, h = 2 h*, A immobile): mean/(V/kappa) is 0.966 +- 0.015 corrected and 1.458 +- 0.023 uncorrected, theory 1 and 1.5. With A in a corner voxel the zero-flux box gives 1.81 / 2.32, so the test puts A at the centre. Tolerance: 5 % plus 3 SEM, and the two must differ by more than 5 % plus 3 combined SEM.
- `docs/ssa.md`: critical spacing with kappa (2 A, fixed reactants), mesh correction, its rules, validation, and the rerun table with both models. README `[ssa]` block and voxel bullets list the key.

## Blockers / Challenges

The example's `radkin/` include prefix refers to the bundled data directory, so only the seed CSV path was made absolute in the scratch copies.

## Commits

- 0a6d17f test: association-time validation of the mesh correction; docs rerun (ticket 03)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
