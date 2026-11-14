# Ticket 06: voxel-rdme-exact

**Acceptance Criteria:**
- [ ] `radkin voxel` with `solver = "ssa"` (config or `--solver`) runs the RDME: per-voxel reactions with the voxel volume, hop events of propensity `D/h²` per molecule and face, zero-flux boundaries, next-subvolume method
- [ ] Initial background: a Poisson draw per voxel and species with mean `c·N_A·V`; Geant4 seeds are injected as exact integer counts at their times
- [ ] Hops and injections conserve molecule numbers exactly (test)
- [ ] Ensemble of N trajectories (`[ssa] trajectories`, `seed`), reduced in trajectory order; outputs keep the existing layout with the ensemble mean, plus `_std` snapshot files and `N_<sp>_std` columns in `yields.csv`
- [ ] `--resume` and checkpoints are rejected with a message under `solver = "ssa"`; `max_coarsen_levels > 0` is rejected until ticket 7
- [ ] A warning is logged when the voxel spacing is below the reaction-validity bound
- [ ] Fixed-seed test on a large-count grid: ensemble-mean species totals agree with the Strang-split LSODA voxel run within 4 SEM plus the documented bias; diffusion alone reproduces the LSODA diffusion of a point seed within 4 SEM
- [ ] ADR 0003 and `docs/ssa.md` describe the method; README documents the voxel usage
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/mod.rs`
- `src/ssa/rdme.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/voxel_output.rs`
- `src/bin/radkin.rs`
- `tests/ssa_rdme.rs`
- `docs/ssa.md`
- `README.md`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including `ssa_rdme`.

**Notes:**

Seed amounts are stored as concentrations (`LoadedSeed::amount`, µmol/L); convert back to integer counts with the fine voxel volume. Species without a diffusion coefficient do not hop.
