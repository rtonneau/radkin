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

The beam comes from the config (`[[beam.<type>]]`, with its `LET`). Logging: `RUST_LOG=info`.

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

## Phase 2 (planned)

Voxelised 3D geometry with species diffusion coefficients, integrated as one ODE system.
rsoda's dense or banded Jacobians do not scale to `n_voxels × n_species` unknowns, so the
Phase-2 design will choose between operator splitting (rsoda per voxel + a diffusion step)
and a sparse/Krylov solver.

## License

LGPL-3.0-or-later, as radiopyo. Bundled data: see [data/README.md](data/README.md).
