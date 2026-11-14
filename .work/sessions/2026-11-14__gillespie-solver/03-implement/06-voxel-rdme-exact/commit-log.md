# Ticket 06: voxel-rdme-exact

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: all `ssa` tests pass, including the 6 unit tests of `ssa::rdme` and the 4 tests of `tests/ssa_rdme.rs`.

`cargo test --release --test ssa_rdme -- --nocapture`:

```
diffusion: worst voxel difference 2.39 SEM
reactions: worst total difference 1.49 SEM, 12 checks above 500 molecules, 9468 reactions and 61929 hops per trajectory
test result: ok. 4 passed; 0 failed; finished in 1.66s
```

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 187 passed, 0 failed.

## Review Notes

Base: 7d218023c930c317ef627790ce41b7199fff7fec. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- `src/ssa/rdme.rs`: next-subvolume method with an indexed binary heap; per-voxel reactions with the voxel-volume model; hops at `D/h^2` per molecule and face, boundary voxels having fewer faces. A hop updates the reaction rate of its two voxels from the reactions that read the hopping species; a full recomputation follows every reaction event, and a rate that drifted above an exact zero produces a null event that resets it.
- Initial solution: Poisson per voxel (`fill_poisson`); seeds as integer counts (`seed_molecules` converts the binned concentrations back and refuses non-integers).
- Conservation: unit test (hops only) and integration test (total exactly 4000 with zero spread at every record).
- Ensemble: `run_rdme`, trajectory i on stream i, reduced in trajectory order; `yields.csv` and the snapshots are byte-identical between 1 and 4 threads (CLI test).
- Outputs: deterministic layout with the mean, `_std` snapshot files, `N_<sp>_std` and `G_<sp>_std` columns; G is the mean of the per-trajectory net count over the deposited energy.
- `--resume` rejected; no checkpoint written; `max_coarsen_levels > 0` and `tau_leap = true` rejected for voxel runs until tickets 7 and 8.
- Small-voxel warning: `RdmeSetup::validity_bound` (unit-tested), logged by `run_rdme`.
- Validation: diffusion field within 2.4 SEM of the deterministic run, reaction totals within 1.5 SEM on a large-count grid.
- Docs: `docs/ssa.md` (RDME section), README (voxel usage); ADR 0003 was written at the grill.

Two things for the user (also in `docs/ssa.md`):
1. On `data/voxel_example.toml` (6 nm voxels) the RDME and the deterministic run differ by up to 30 % on G-values (e_aq 1.13 vs 0.84, H2O2 1.77 vs 2.02). This is the low-count regime the Idea was written for, not a defect: with about one molecule per voxel the rate equations make a molecule react with itself.
2. The constant 0.2527 of the validity bound `h* = 0.2527 k / (N_A D)` and its reference (Hellander, Hellander and Petzold 2012) are quoted from memory and were not checked against the paper. The bundled chemistry triggers the warning at 6 nm because reactions written with H2O as a reactant get an effective constant `k [H2O]` (75 nm for `e_aq + O_minus + H2O`).

## Blockers / Challenges

The shell heredocs of this environment collapse backslash sequences, which broke two string literals during an edit; they were repaired with the editor tool. No blocker.

## Commits

- 3f1be8e feat: voxel RDME with hop events, Poisson background and Geant4 seeds (ticket 06)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
