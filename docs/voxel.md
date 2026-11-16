# Voxel reaction-diffusion

`radkin voxel CONFIG --out DIR` runs the chemistry of `CONFIG` (with its includes) on a
uniform Cartesian voxel grid seeded by Geant4-DNA molecules. Each splitting step is
Strang: reaction(dt/2), diffusion(dt), reaction(dt/2) (ADR 0002 in `.work/adr/`).

- The reaction half-step is an independent rsoda run per voxel, parallelised with rayon
  (cargo feature `parallel`, on by default; the results are bit-identical to serial).
- Diffusion uses Crank-Nicolson ADI with zero-flux boundaries. After each injection,
  backward-Euler substeps (Rannacher start) keep the concentrations non-negative.
- dt is adaptive by step doubling and lands exactly on seed and snapshot times.
- G-values in `yields.csv` are net of the background chemistry: the seedless uniform
  solution is integrated alongside, so a background that drifts (e.g. a pH slightly off
  equilibrium) does not show up as yield.
- Background voxels, which still hold only the uniform initial solution (within
  `1e-3 * solver_atol + 1e-12 * |c|`), share one rsoda run per half-step, so untouched
  voxels cost almost nothing. Diffusion then acts on the deviation from that background.
  Set `background_voxels = false` in `[voxel.run]` to integrate every voxel.
- Global coarsening (off by default) merges every 2x2x2 block of voxels into one once the whole
  field is smooth, and refines back to the fine grid at the next injection. It is meant for runs
  that continue until the system has homogenised, where it cuts the cost of the late steps. It
  does nothing while a spur is still localized (see the `[voxel]` keys below).

```toml
[includes]
toml = ["radkin/CONFIG_PAPER_01.toml"]   # chemistry

[voxel]
clip = false                 # true: drop out-of-grid molecules (lost fraction is logged)
# max_coarsen_levels = 0     # halvings allowed (0 = off); each grid dim must be divisible by 2^levels
# coarsen_tol_factor = 0.1   # > 0: coarsen when the spread is below factor * (atol + rtol * |mean|)

[voxel.grid]
dims = [8, 8, 12]
spacing_nm = 6.0             # or [dx, dy, dz]
origin_nm = [-24, -24, -36]

[voxel.diffusion]            # m^2/s; defaults are Geant4-DNA option3's (HO3: O3's)
O2 = 2.4e-9                  # (e_aq, OH_r, ..., O_atom); other species do not diffuse

[voxel.geant4]
aliases = { "O2(B)" = "O2" } # Geant4 name -> radkin label (option3 species built in)
ignore = ["H2O^0"]

[[voxel.seed]]               # one entry per injection
file = "track.csv"           # species,x,y,z[,event_id]; positions in nm
time = 1e-12                 # s
events = [0]                 # optional event_id filter
offset_nm = [0, 0, 0]
energy_ev = 360.0            # optional: deposited energy, enables G(t)

[voxel.run]
t_end = 1e-7
snapshots = [1e-12, 1e-9, 1e-7]
rtol = 1e-2                  # splitting error control (atol in µmol/L)
atol = 1e-3
# background_voxels = true   # share one integration for untouched voxels
# checkpoint_every = 2       # checkpoint every 2 snapshot times (default: end only)
```

## Coarsening keys (`[voxel]`)

- `max_coarsen_levels` (integer >= 0, default 0 = off): how many times the grid may be halved
  (dims / 2, spacing x 2 per level). Every dimension of `[voxel.grid] dims` must be divisible
  by `2^max_coarsen_levels`; a violation is a configuration error naming the axis.
- `coarsen_tol_factor` (> 0, default 0.1): a level is dropped when, for every species and every
  2x2x2 block, the spread (max - min) over the block and its face neighbours is at most
  `coarsen_tol_factor * (atol + rtol * |mean|)`, with `atol` and `rtol` of `[voxel.run]`. Smaller
  is more conservative. On the example case at `t_end = 1e-5` the domain totals stay within
  1e-6 (relative) of the fine run up to a factor of about 3; see `docs/perf/voxel.md`.

Merging conserves the amount of every species, and diffusion and reaction then run on the
coarse voxels. Snapshots and checkpoints are always written on the fine grid (the coarse state
is copied to its fine voxels) and the `yields.csv` sums are volume-weighted, so outputs have the
same shape with or without coarsening. A resumed run recomputes its level from the checkpoint. The coarsened result agrees
with the fine run within the merge tolerance, not bit for bit.

The background species come from the chemistry's initial concentrations, uniform over the
grid. Each molecule adds `1/(N_A V)` to its voxel. The molecule CSV is written on the
Geant4 side, e.g. from a chem `TimeStepAction` with `G4MolecularConfiguration::GetName()`
(`e_aq^-1`, `°OH^0`, `H3O^1`, ...) or the user IDs (`e_aq`, `°OH`, `H3Op`, `OHm`, ...).

## Outputs

| file | content |
|---|---|
| `yields.csv` | `t`, `N_<sp>` (molecules in the domain) at every accepted step; with `energy_ev`, `G_<sp>` = (N - N(t_start)) / (E / 100 eV) |
| `snapshots/t<k>_<sp>.npy` | concentration field [mol/L], C order, shape `(nz, ny, nx)` |
| `snapshots.json` | grid, species, snapshot times and files, step counts |

## Checkpoints

Every voxel run writes `DIR/checkpoint.npz` at `t_end`, and every
`checkpoint_every` snapshot times if set (each write replaces the file atomically).
`radkin voxel CONFIG --resume OLD/checkpoint.npz --out NEW` starts from it instead of
`t_start`, for three uses:

- recovery: rerun the same configuration after a crash or a kill;
- extension: raise `t_end` and continue a finished run;
- new conditions: change the beam, solver settings or seeds. Seeds already injected are
  skipped; new seeds must come at or after the checkpoint time.

The grid (dims, spacing, origin) and the species order must match the checkpoint. Other
differences are logged as warnings. The resumed run restarts its step control from the
saved dt, so it agrees with an uninterrupted run within the splitting tolerance, not bit
for bit. It writes only the part after the checkpoint into `NEW` (`N(t_start)` of its
G-values is taken at the checkpoint time); joining the outputs is left to the user. The
checkpoint is a NumPy archive: `np.load("checkpoint.npz")` gives `state` (µmol/L, shape
`(n_species, nz, ny, nx)`), `time` (`[t, dt]`), `grid`, `seeds`, `background` and
`meta.json` (species, dims, settings summary).

`python -I tools/compare_geant4.py DIR [--geant4 g4.csv]` plots G(t), optionally against
Geant4-DNA chemistry yields (a CSV with columns `t,species,G`).

## Mean-field caveat

Each voxel is well mixed, so one molecule in a (5 nm)^3 voxel is
about 13 mM, and second-order terms such as `e_aq + e_aq` act even on a lone molecule.
Voxels larger than the spur size (a few nm) dilute a spur over the voxel and
underestimate recombination: on a 100 keV electron track (Geant4-DNA option3 seeds at
1 ps), G(OH) at 10 µs is 5.2 with 1 µm voxels, 4.3 with 250 nm and 4.2 with 40 nm on a
track segment, against 2.1–3.3 for Geant4-DNA IRT. Voxels much smaller than the mean distance between molecules make
the continuum description meaningless. Check that the results are stable as the voxel
size changes, and compare against Geant4-DNA step-by-step chemistry before trusting
absolute G-values. The example (`data/voxel_example.toml`, a synthetic track) shows the
workflow; its numbers are illustrative only.
