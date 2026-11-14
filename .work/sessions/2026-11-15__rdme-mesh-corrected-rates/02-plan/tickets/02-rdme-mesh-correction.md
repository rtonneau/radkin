# Ticket 02: rdme-mesh-correction

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] A helper gives the pair constant `kappa` [m^3/s] and summed `D` of a mass-action reaction with two dynamic slots: `kappa = c * N_A V / (1e3 N_A)`, times 2 for a single species with multiplicity 2; `D` summed with multiplicity.
- [ ] `validity_bound` uses it: `h* = 0.2527 kappa / D` (the `2 A` value doubles; its unit test updated and a `2 A` case added).
- [ ] `RdmeSetup::new` takes the resolved correction; with Hellander on a cubic grid, each level with spacing `h` scales `c` by `1/(1 - h*/h)` when `h > h*` and leaves it unchanged otherwise; Michaelis, order != 2 and `D = 0` reactions are untouched; fixed reactants enter through `c` (k*[fixed]).
- [ ] `run_rdme` resolves the setting: unset or explicit Hellander on a cubic grid → corrected; unset on an anisotropic grid → warning, uncorrected; explicit Hellander on an anisotropic grid → config error; `"none"` → uncorrected. One warning per reaction left uncorrected at `h <= h*` on the fine grid (h, h*, reaction string), replacing the existing below-`h*` warning when the correction is on; the existing warning stays for `"none"`.
- [ ] Unit tests: the factor on levels 0 and 1 with their own `h`; the pole fallback; the `2 A` kappa; a fixed-reactant reaction; `"none"` equals the old model.
- [ ] `cargo test --release` passes (existing RDME integration tests with the new default).

**Files to Touch:**
- `src/ssa/rdme.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/rdme_coarsen.rs` (test setup call)
- `src/ssa/rdme_tau.rs` (test setup call)

**Verification Step:**

Run:
```bash
cargo clippy --release --all-targets -- -D warnings && cargo test --release ssa
```

Expected:
No warnings; all `ssa` tests pass, including the new mesh-correction unit tests.

**Notes:**

`h* = C_3 kappa / (6 D)` with C_3 = 1.5164 (0.2527) comes from Hellander & Hellander 2020 eq. 6-8 with `1/k = 1/k_r + 1/(4 pi sigma D)`: `k_meso h^3 = k / (1 - h*/h)`. Cubic means the three spacings equal within a relative 1e-12. The per-level `h` is `grid.spacing[0]` of that level.
