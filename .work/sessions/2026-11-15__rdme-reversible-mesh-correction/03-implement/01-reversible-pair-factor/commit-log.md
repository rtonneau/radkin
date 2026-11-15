# Ticket 01: reversible-pair-factor

**Status:** ✅ Done

## Local Test Result

Verification Step, run from the repo root (log: `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/evidence/verification-ssa.log`):

```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa -- --nocapture
```

Exit 0: fmt clean, clippy clean, every `ssa` test passes (lib 67, ssa_rdme 10, ssa_rdme_multigrid 3, the other ssa_* files unchanged). Full `cargo test --release`: 235 passed, 0 failed, 1 ignored (the existing `ssa_aerated_track_measurement`) (`evidence/verification-full.log`).

New seedless pH 7 background tests (CONFIG_PAPER_01, 32^3 voxels, 10 us, mean change of the ions against 20 000 unit-cell SSA trajectories of the box volume; per size |z| <= 4 and the pooled `sum z / sqrt(n)` <= 4). H3O_plus shown, OH_minus is identical (only the water pair runs):

| test | 5 nm | 10 nm | 20 nm | pooled |
|---|---|---|---|---|
| single grid, after | +0.67 | +0.19 | -0.78 | +0.05 |
| single grid, before (base 4ea5947) | -3.59 | -3.57 | -3.01 | **-5.87 (fails)** |
| multi-grid, OH- level 0, after | +0.93 | -0.06 | +0.07 | |
| multi-grid, OH- level 1, after | -0.86 | -0.10 | -0.64 | -0.27 (all six) |
| multi-grid, OH- level 0, before | -3.02 | -3.97 | -2.49 | |
| multi-grid, OH- level 1, before | -2.58 | -2.19 | -0.92 | **-6.19 (fails)** |

"Before" ran the new test files unchanged on a `git archive` of the base commit in the scratch dir (they use only the pre-existing API): `evidence/before-fix-background.log`, both tests fail on the pooled check. Wall time on this machine (16 threads): about 10 s single grid, 14 s multi-grid.

Setup logging checked with `RUST_LOG=info radkin voxel ... --solver ssa` on a 4^3, 5 nm CONFIG_PAPER_01 box plus a made-up pair bimolecular in both directions: one info line per corrected pair, one warning for the K_eq shift (`evidence/setup-log.txt`).

## Review Notes

Checked against the Acceptance Criteria:

- AC1 (seedless box test, fails before, passes after, single grid and multi-grid): `ssa_rdme_background_keeps_the_water_equilibrium` (tests/ssa_rdme.rs) and `ssa_multigrid_background_keeps_the_water_equilibrium` (tests/ssa_rdme_multigrid.rs), shared helper `tests/common/mod.rs`. Results in Local Test Result. The drift over 10 us is only 1.5-6 % of the ions (relaxation limited by the encounter kinetics, not the full 1/(1-h*/h)), so one size alone cannot reach 4 SEM at an affordable cost; each size is checked at 4 SEM and the pooled z of the sizes is checked at 4 SEM too, which is what fails before the fix. The tests use the low-level API (`RdmeSetup`/`Rdme`/`Multigrid` + `fill_poisson`) rather than `run_rdme`, about 3x cheaper. O2 is kept at its CONFIG_PAPER_01 concentration but does not diffuse (nothing reacts with it in the box; its hops would dominate). The multi-grid test puts every non-ion species on level 2 (cheap fine voxels) and runs the ions on level 0, then OH- on level 1 (the cross-level case).
- AC2 (pairs readable, unit tests on the radiopyo network): `SimMatrices::reversible_pairs` `(forward, backward)` in `k` order, and `SsaModel::pairs` + `SsaModel::is_bimolecular`. Tests: `matrices::tests::reversible_pairs_are_the_declared_pairs` (the 9 declared CONFIG_PAPER_01 pairs; the one-way `e_aq + H2O -> H_r + OH_minus` is not a pair) and `ssa::model::tests::reversible_pairs_and_their_bimolecular_direction` (the bimolecular direction of each).
- AC3 (shared factor in `mesh_correct`, per level, multi-grid levels; both-bimolecular pairs separate + log): `RdmeSetup` classifies pairs (`ReversiblePair::Shared { bimolecular, follower }` / `BothBimolecular`) and `mesh_correct` divides the follower by the `1 - h*/h` of its bimolecular direction on every level. In the multi-grid the two directions can fire on different levels (forward of water on the finest product level, backward on the coarsest reactant level), so each `Layout` now holds per-level propensity models where the follower carries the factor of the level where its bimolecular direction fires; the propensity call sites in `rdme_multigrid.rs` and `rdme_multigrid_tau.rs` read `lay.models[l]`. Unit tests: `reversible_pairs_share_the_factor_of_their_bimolecular_direction`, `reversible_pair_follower_stays_plain_at_the_pole` (rdme.rs), `a_reversible_pair_takes_the_factor_of_its_bimolecular_level` (rdme_multigrid.rs). Unaffected reactions keep bit-identical constants (division by the same divisor as before, layout models are clones).
- AC4 (logged once at setup): `log_reversible_pairs` in `voxel_sim::run_rdme`, called once after `RdmeSetup::new` when the correction is on: info per shared pair with its fine-grid factor, warning with the K_eq shift for both-bimolecular pairs.
- AC5: fmt, clippy -D warnings, full `cargo test --release` pass.

Scope notes: `src/ssa/rdme_multigrid_tau.rs` (3 propensity call sites) and `src/ssa/mod.rs` (re-export `ReversiblePair`) are outside the listed files but needed; `tests/common/mod.rs` is new (shared by the two test files). Every hellander voxel SSA result with a reversible pair changes (expected, ADR 0006); the existing ssa tests still pass. The kernel and the deterministic/unit-cell paths are untouched. No review rounds yet.

## Blockers / Challenges

- Statistical power vs test time: the 10 us drift is small (see Review Notes), and per-trajectory cost is dominated by the two full voxel rebuilds of `Rdme::new` + `fill_poisson` (~60 propensities per voxel), not by events. The pooled check keeps the two tests at about 24 s together on 16 threads; on the 4-core CI runners expect roughly 4x that. Lowering it further would need fewer trajectories (weaker pre-fix failure) or a rebuild that skips empty voxels (product change, out of scope).
- Python edits briefly wrote CRLF into a few files; normalised back to LF (`git ls-files --eol` clean).

### Parent review and acceptance

- Base 4ea59473fcb126d0b9b51f3c03a5331bc659a87f; pre-existing changes: gps session state only. Recovery/correction passes used: 1.
- Independent review (fresh reviewer, base 4ea5947): no P0-P2. P3-1 `src/ssa/voxel_sim.rs` log_reversible_pairs logged "corrected" at/below the pole: fixed by the parent (guard `pole_distance(0, bimolecular) < 1.0`). P3-2 `tests/ssa_rdme_multigrid.rs` correlated z values across OH- levels: fixed by the parent (seed `31 + oh_level`, plus a comment on the shared unit-cell reference).
- Re-review of the corrected diff: no findings. Its gap (fails-before not re-shown with the new seeds) is closed: the current test on a base archive fails at pooled -5.59 SEM (`evidence/01-before-fix-new-seeds.log`). The test comment was updated to these numbers.
- Verification after the fixes: fmt, clippy and `cargo test --release ssa` pass (`evidence/01-parent-fix.log`).
- Accepted: every Acceptance Criterion is met and verified, and no actionable findings remain. Requested model opus/high for the implementer and the reviewer.
- Remaining gap: CI runtime on 4-core runners and Ubuntu results are not checked locally.

## Commits

- 084cfc6 fix: share the mesh-correction factor across reversible pairs (ticket 01)

## Time Spent

44m (ticket-start.js to ticket-complete.js)
