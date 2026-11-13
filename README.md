# radkin

Homogeneous chemical kinetics of water radiolysis and radiobiology under conventional and
ultra-high-dose-rate irradiation, in Rust. radkin is a port of
[radiopyo](https://github.com/rtonneau/radiopyo) (Python) and reproduces its results; the
ODE system is integrated with [rsoda](https://github.com/rtonneau/rsoda), a Rust port of
LSODA, with radiopyo's analytic Jacobian.

Same model as radiopyo: G-values (Boscolo 2020 LET fit, optional dose-rate factor) feed
rate equations for the species concentrations; acid/base, equilibrium, k-, Arrhenius,
constant and Michaelis–Menten reactions written in TOML; constant, single-pulse, pulsed and
macro-pulsed beams; clipped concentrations, smooth depletion switch and overshoot
relaxation (radiopyo ADRs 0001–0008). Configuration files of radiopyo work unchanged
(`radiopyo/` include prefix included).

## Build

```bash
cargo build --release      # target/release/radkin
cargo test --release       # unit tests + golden comparisons with radiopyo
```

Rust 2024 edition (stable toolchain). No C dependency.

## Command line

```bash
radkin run data/config_Labarbe.toml --t-span 1e-7 10 --out out/labarbe
radkin run data/tonneau_2025.toml --mode single-pulse --t-span 0 500 --out out/tonneau
radkin run my_pulsed.toml --mode chunked --t-span 0 1 \
    --max-step-on 1e-7 --max-step-off 1e-3 --points-on 10 --points-off 50 \
    --production H2O2,O2 --consumption O2 --out out/pulsed
```

| Option | radiopyo equivalent |
|---|---|
| `--mode run` (default) | `UnitCell.run(t_span, ...)` |
| `--mode chunked` | `prepare_chunked_run(...).run()` |
| `--mode single-pulse` | `prepare_single_pulse(...).run()` |
| `--t-eval FILE` | `t_eval=` (times separated by blanks; `run` mode) |
| `--rtol --atol --max-step --first-step --no-jac` | `rtol`, `atol`, `max_step`, `first_step`, `use_jac=False` |
| `--sparse-jacobian` | none: sparse Jacobian and sparse LU (rsoda feature `sparse`) instead of dense |
| `--max-step-on/off/trans/final`, `--trans-time` | `max_step_size_on/off/trans/final`, `trans_time` |
| `--points-on/--points-off N` | `points_per_chunk_on/off` |
| `--beam LABEL`, `--let X`, `--temperature C` | `set_beam(label, LET=)`, `change_LET`, `set_temperature` |

`--sparse-jacobian` (library: `RunOptions::jacobian = JacobianMode::Sparse`)
builds the Jacobian's sparsity pattern from the reaction tables and solves the
stiff steps with a sparse LU. Results match the dense run within the golden
tolerances. For the single cell (15 to 22 species) it is not faster; it is
meant for larger systems. `cargo run --release --example bench -- CONFIG
--t-span 0 1 --jacobian both` times both.

The beam comes from the config (`[[beam.<type>]]`, with its `LET`).
`data/tonneau_2025.toml` is the chemistry of Tonneau et al. 2025 (PMB 70 235021) with the
inputs of the paper's original runs (variable G-values, LET 0.56 keV/µm, 30 Gy single
pulse at 1e7 Gy/s). Logging: `RUST_LOG=info`.

## Output files

Each run writes into `--out`:

| File | radiopyo | Units |
|---|---|---|
| `concentrations.csv` | `res.to_pandas()` | µmol/L; `t` then dynamic species, then fixed species |
| `reaction_rates.csv` | `res.reconstruct_reactions()` | mol/L/s; columns `KR_i`/`CR_i` (pairs first), `G_<sp>`, `EZ_i` |
| `species_integrals.csv` | `res.integrate_species()` | µmol/L·s |
| `production_<sp>.csv` | `res.integrate_production_of(sp)` | µmol/L, sorted descending |
| `consumption_<sp>.csv` | `res.integrate_consumption_of(sp)` | µmol/L, sorted descending |
| `run.json` | — | inputs, G-values, bursts, short ON spans, solver counters, `label -> str(reaction)` |

`pd.read_csv("concentrations.csv", index_col=0)` gives the same DataFrame as
`res.to_pandas()`.

## Library

```rust
use radkin::unit_cell::{RunOptions, UnitCell};

let mut cell = UnitCell::from_toml("data/config_Labarbe.toml")?;
let res = cell.run((1e-7, 10.0), &RunOptions::default(), None)?;
let (labels, rows) = res.concentrations();
```

The kernel (`radkin::kernel::{derive, jacobian, reaction_rates}`) works on one cell's state
slice with caller-provided buffers: it is meant to be called per voxel by the Phase-2
spatial driver.

## Comparison with radiopyo

`tools/` holds Python scripts run in an environment where radiopyo is importable:

- `export_golden.py`: exports radiopyo's golden data to `tests/golden/*.json` (committed).
- `compare_radiopyo.py`: runs both codes on the same inputs (a wrapper TOML including your
  config and adding the beam) and compares species, reaction rates and integrals.
- `acceptance.ps1`: the Phase-1 matrix — pure water, Labarbe and the Tonneau 2025 paper
  config, × constant, single-pulse, pulsed and macro-pulsed beams, × run / chunked
  (/ single-pulse) modes. All 27 cases pass at `rtol = 1e-6`; the observed differences are
  ~1e-8 relative.

```powershell
pwsh tools/acceptance.ps1 -Python <python with radiopyo> -Radiopyo <radiopyo checkout>
```

Results agree to solver tolerance, not bit for bit: scipy's `solve_ivp` wraps ODEPACK's
LSODA while rsoda is a direct port of the thread-safe Fortran fork.

## Voxel reaction-diffusion (Phase 2)

`radkin voxel CONFIG --out DIR` runs the chemistry of `CONFIG` (with its includes) on a
uniform Cartesian voxel grid seeded by Geant4-DNA molecules. Each splitting step is
Strang: reaction(dt/2), diffusion(dt), reaction(dt/2) (ADR 0002 in `.work/adr/`).

- The reaction half-step is an independent rsoda run per voxel, parallelised with rayon
  (cargo feature `parallel`, on by default; the results are bit-identical to serial).
- Diffusion uses Crank-Nicolson ADI with zero-flux boundaries. After each injection,
  backward-Euler substeps (Rannacher start) keep the concentrations non-negative.
- dt is adaptive by step doubling and lands exactly on seed and snapshot times.
- Background voxels, which still hold only the uniform initial solution (within
  `1e-3 * solver_atol + 1e-12 * |c|`), share one rsoda run per half-step, so untouched
  voxels cost almost nothing. Diffusion then acts on the deviation from that background.
  Set `background_voxels = false` in `[voxel.run]` to integrate every voxel.

```toml
[includes]
toml = ["radkin/CONFIG_PAPER_01.toml"]   # chemistry

[voxel]
clip = false                 # true: drop out-of-grid molecules (lost fraction is logged)

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
```

The background species come from the chemistry's initial concentrations, uniform over the
grid. Each molecule adds `1/(N_A V)` to its voxel. The molecule CSV is written on the
Geant4 side, e.g. from a chem `TimeStepAction` with `G4MolecularConfiguration::GetName()`
(`e_aq^-1`, `°OH^0`, `H3O^1`, ...) or the user IDs (`e_aq`, `°OH`, `H3Op`, `OHm`, ...).

Outputs:

| file | content |
|---|---|
| `yields.csv` | `t`, `N_<sp>` (molecules in the domain) at every accepted step; with `energy_ev`, `G_<sp>` = (N - N(t_start)) / (E / 100 eV) |
| `snapshots/t<k>_<sp>.npy` | concentration field [mol/L], C order, shape `(nz, ny, nx)` |
| `snapshots.json` | grid, species, snapshot times and files, step counts |

`python -I tools/compare_geant4.py DIR [--geant4 g4.csv]` plots G(t), optionally against
Geant4-DNA chemistry yields (a CSV with columns `t,species,G`).

**Mean-field caveat.** Each voxel is well mixed, so one molecule in a (5 nm)^3 voxel is
about 13 mM, and second-order terms such as `e_aq + e_aq` act even on a lone molecule.
Voxels larger than the spur size (a few nm) mix a spur too early and overestimate
recombination. Voxels much smaller than the mean distance between molecules make
the continuum description meaningless. Check that the results are stable as the voxel
size changes, and compare against Geant4-DNA step-by-step chemistry before trusting
absolute G-values. The example (`data/voxel_example.toml`, a synthetic track) shows the
workflow; its numbers are illustrative only.

## License

LGPL-3.0-or-later, as radiopyo. Bundled data: see [data/README.md](data/README.md).
