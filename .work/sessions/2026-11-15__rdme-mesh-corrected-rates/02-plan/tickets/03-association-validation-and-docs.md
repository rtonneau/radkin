# Ticket 03: association-validation-and-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] An integration test `ssa_rdme_mesh_correction_matches_association_time` (in `tests/ssa_rdme.rs`): one A and one B, `A + B -> C`, uniform random start, cubic grid with `h` a small multiple of `h*`; the mean first-reaction time with the correction matches `V/k` within a tolerance derived from its SEM, and without it is longer by more than that tolerance (Hellander 2012 eq. 1).
- [ ] `docs/ssa.md` documents `rate_correction` (formula, default, `h <= h*`, fixed reactants, `2 A`, anisotropic and unit-cell rules) and fixes the critical-spacing paragraph (`2 A` pair constant).
- [ ] The small-voxel table of `docs/ssa.md` is rerun on `data/voxel_example.toml` (64 trajectories, seed 1) with and without the correction, with the deterministic column kept.
- [ ] `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` are clean.

**Files to Touch:**
- `tests/ssa_rdme.rs`
- `docs/ssa.md`
- `README.md` (only if it lists `[ssa]` keys)

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks clean; the association-time test passes.

**Notes:**

Run outputs go to `.scratch/tests/2026-11-15__rdme-mesh-corrected-rates`. The mean association time for one pair in volume `V_tot` is `V_tot/kappa` when the macroscopic model holds. Use `D_A = 0` (B diffuses) to keep it cheap, and a reflective-boundary `tau_D` check if the measured gap needs explaining.
