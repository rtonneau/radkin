---
bump: minor
floor: minor
---
- Stochastic solver: `solver = "ssa"` in the run file, or `--solver ssa` on `radkin run` and `radkin voxel`, runs the chemistry on molecule counts with Gillespie's algorithm instead of the rate equations. It is configured by a new `[ssa]` table (`volume`, `trajectories`, `seed`, `tau_leap`, `tau_epsilon`).
- `radkin run --solver ssa` writes the ensemble mean in the usual `concentrations.csv` layout, followed by one `<species>_std` column per species. A seed gives the same output for any number of threads.
- `radkin voxel --solver ssa` runs the voxel grid stochastically: reactions in each voxel and hops of molecules between voxels. Outputs keep their layout and gain `_std` snapshot files and `N_<sp>_std`, `G_<sp>_std` columns in `yields.csv`. `max_coarsen_levels` is supported; checkpoints and `--resume` are not.
- `tau_leap = true` switches to tau-leaping, an approximation that is much faster on large molecule counts without fast equilibria; `tau_epsilon` sets its accuracy.
- `docs/ssa.md` describes the solver and its validation against LSODA; `examples/ssa_bench.rs` compares the exact and tau-leaping runs.
