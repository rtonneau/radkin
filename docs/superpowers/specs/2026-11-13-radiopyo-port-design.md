# radkin — Rust port of radiopyo (Phase 1)

Date: 2026-11-13
Status: draft, awaiting approval

## Intent

`radkin` is a Rust translation of radiopyo (`D:\DEV\PYTHON\radiopyo`,
GitHub `rtonneau/radiopyo`, LGPL-3.0): homogeneous chemical kinetics of
water radiolysis and radiobiology under CONV / UHDR irradiation. The ODE
system is integrated with the local Rust LSODA port `rsoda`
(`D:\DEV\RUST\rsoda`, public GitHub `rtonneau/rsoda`) instead of
`scipy.integrate.solve_ivp`.

Why: Phase 2 extends the model to a voxelised 3D geometry where species
diffuse between voxels via diffusion coefficients (a DSMC-like mesoscopic
picture, no individual trajectories), still integrated as an ODE system.
That needs a fast, compiled core. Phase 1 is the trustworthy base: a port
that reproduces radiopyo's numbers.

## Goals (Phase 1)

1. Reproduce radiopyo's results on the bundled configurations and every
   beam type, within solver tolerance (see "Acceptance").
2. Produce data files that compare directly with radiopyo's pandas
   outputs (layout A, below).
3. Keep the chemistry kernel (rates, right-hand side, Jacobian) independent
   of how many "cells" the state vector holds, so Phase 2 can reuse it per
   voxel.

## Non-goals (Phase 1)

- Voxel geometry and diffusion (Phase 2, separate design).
- HDF5 persistence / compatibility with radiopyo's `H5Manager`,
  `SimCollection`, `SimIdentifier`.
- Python bindings.
- pandas-style analysis helpers beyond the files listed under "Outputs"
  (unit converters `O2Convert` etc., `ResCollection`, equilibrium-yield
  helpers, `species_balance_evolution` — the latter may be added if cheap).
- Other solve_ivp methods (BDF, Radau, RK45...). Only LSODA.
- Logging-config parity. Use the `log` crate (or `tracing`) with warnings
  where radiopyo logs warnings.

## Anti-goals

- Silently "fixing" radiopyo's chemistry. Every modelling convention of
  radiopyo is reproduced as-is, including its ADRs (see "Semantics to
  preserve"). Divergences are bugs, unless documented and approved.
- Returning a truncated result when the solver fails (radiopyo raises
  `SolverFailedError`; radkin returns an error).

## Constraints

- Rust edition 2024 (same as rsoda), stable toolchain.
- rsoda as a git dependency on `https://github.com/rtonneau/rsoda`
  (pinned rev), so CI builds without the local checkout. A local
  `[patch]` override to `../rsoda` is allowed for development (documented,
  not committed as the default).
- rsoda API: `Lsoda::new(system, options, t0, y0)`, `OdeSystem { rhs, jac,
  has_jac }`, dense Jacobian column-major `pd[i + j*neq]`,
  `OptionsBuilder` with `tolerance`, `jacobian(JacobianKind::DenseUser)`,
  `h0`, `hmax`, `mxstep`; `integrate`, `step`, `integrate_tcrit`,
  `step_tcrit`; interpolation module for dense output.
- radiopyo is LGPL-3.0, so radkin is LGPL-3.0. The bundled data files
  (`advanced_water_radiolysis.toml`, `reactions_biology.toml`,
  `config_Labarbe.toml`, `gvalues_from_Boscolo_2020.csv`) are copied into
  `radkin/data/` with attribution.
- Private GitHub repo `rtonneau/radkin`, project root
  `D:\DEV\RUST\radkin`.
- No C library dependencies (keeps Windows + CI simple).

## Approaches considered (output format)

| | Approach | Trade-off |
|---|---|---|
| **A (chosen)** | CSV files mirroring radiopyo's DataFrames + a Python compare script | No heavy deps; one-to-one with `to_pandas()` / `reconstruct_reactions()`; `pd.read_csv(index_col=0)` gives the same DataFrame |
| B | HDF5 readable by radiopyo `SimCollection` | Needs HDF5 C library (painful on Windows/CI); must replicate `SimIdentifier` naming exactly |
| C | NPZ | Easy for numpy, loses readable labels, needs an extra crate |

The user picked A.

## Scope of the port (what must exist)

### 1. Configuration (radiopyo `parser/`)

- TOML loading with `[includes]` resolved recursively, exactly as
  radiopyo: an entry starting with exactly `radiopyo/` maps to the bundled
  data directory (keep this prefix so radiopyo configs work unchanged;
  `radkin/` may be added as an alias); relative paths resolve against the
  including file's folder, then the CWD; absolute paths accepted.
- Merge rules of `parser/utils.py` (`combine_config`, `reaction_key`):
  precedence (includes in list order, then the including file), entry-wise
  merge of `[parameters]`, `[concentrations.fixed]`,
  `[concentrations.initial]`, `[reactions.radiolytic]`; reactions keyed on
  their species on each side, order-insensitive; field-by-field override;
  the `equilibrium` exception (>= 2 of the 3 constants given → fully
  redefined); redefined reactions keep their place, new ones are appended;
  beam rules (the `default` label is never imported from includes).
- Sections: `[parameters]` (`pH`, `force_cst_pH`, `temperature`,
  `variable_GValues`, `species_order`), `[concentrations.fixed]`,
  `[concentrations.initial]`, `[reactions.radiolytic]`,
  `[[reactions.acid_base]]`, `[[reactions.equilibrium]]`,
  `[[reactions.k_reaction]]`, `[[reactions.constant]]`,
  `[[reactions.arrhenius]]`, `[[reactions.enzymatic]]`, beams
  (`constant`, `pulsed`, `singlepulse`, `macro_pulsed`, with every
  alternative parametrisation: `frequency`, `micro_frequency`,
  `pulse_dose`, `n_pulse` / `max_dose`, `start_time`, `micro_start_time`).
- Reaction-string grammar of `parser/grammar/reaction_grammar.lark`
  (k, constant `-> X`, equilibrium `<->`, enzymatic `A -- enzyme >> ...`,
  multiplicity, signed "negative" products, charges `_plus`/`_minus`/`+`/`-`,
  radicals `_r`, `e_aq`, plain-name elements). A hand-written parser is
  fine; no need to replicate Lark. Error messages should name the
  offending string.
- k values: scalar, list (log10 polynomial of 1/T), Arrhenius
  (`A`, `Ea` in kJ/mol), equilibrium constants. Temperature from
  `[parameters].temperature` (degC), the single source of truth.

### 2. Chemistry (radiopyo `chemistry/`)

- Species list in deterministic order of appearance (or `species_order`),
  split into dynamic (solved) and fixed species. Fixed species and initial
  concentrations may be temperature polynomials (`cc_init`).
- Reaction kinds: k-reaction, constant, Arrhenius, acid/base (two
  k-reactions), equilibrium (two k-reactions, third constant deduced),
  Michaelis–Menten enzymatic, radiolytic (G-value production).
- Reaction labels identical to radiopyo `as_label()` (needed for CSV
  column names).
- Global reaction index order of `ReactionsList.__post_init__`: acid/base
  (2 each), equilibrium (2 each), k-reactions, enzymatic, radiolytic.
- pH: initial H3O+/OH- from pH; `force_cst_pH` as in radiopyo.

### 3. Matrices and right-hand side (radiopyo `environment_build.py`,
`ode_functions_solver.py`)

- Equivalent of `SimMatrices` (`KMatrices`, `GVector`, `EzMatrices`,
  `StateLayout`, `OvershootTables`), built once per environment/beam.
- `reaction_rates(t, y_full, mat, dose_rate)` and `derive` / `jacobian`
  with the same math. State in µmol/L for dynamic species, internal molar
  full state with fixed species appended (as `StateLayout.to_full_molar`).
- **Kernel shape for Phase 2:** the rate/RHS/Jacobian functions take one
  cell's state slice and write into caller-provided slices (no
  allocation per call, no global state), so a Phase-2 driver can call them
  per voxel.

### 4. Beams (radiopyo `physics/`)

- `Timer` / `MacroTimer` state machines, `ConstantBeam`,
  `SinglePulseBeam`, `PulsedBeam`, `MacroPulsedBeam`, `BeamCollection`,
  the generic transition walk (`make_matrix`), `at(t).dose_rate()`,
  `MIN_BEAM_ON_TIME` = 1e-6 s → `ShortBeamOn` error at build time.

### 5. G-values (radiopyo `utils/g_values.py`)

- Boscolo 2020 interpolation from the bundled CSV, valid LET 0.13–150
  keV/µm; rewrites e_aq, OH_r, H_r, H2, H2O2 and H3O_plus (= e_aq, ADR
  0002); optional `variable_GValues` dose-rate factor
  `0.97207704 * dose_rate_peak^-0.10972` without range check. `set_beam`
  without LET reuses the previous LET.

### 6. Driver (radiopyo `unit_cell.py`)

- `UnitCell` equivalent with `run(t_span, options)`,
  `run_till_max_dose`, `prepare_chunked_run(...).run()`,
  `prepare_single_pulse(...).run()`, `set_beam`, `set_temperature`.
- `RunOptions`: `atol` (default 1e-12), `rtol` (default 1e-10),
  `max_step`, `first_step` (with radiopyo's reduction rules:
  `>= max_step` → `max_step/10`; `> span` → `span/1000`), `use_jac`
  (default true: analytic dense Jacobian; false → rsoda finite
  differences), `t_eval`, and for chunked/single-pulse runs
  `max_step_size_on/off/trans`, `points_per_chunk_on/off` (point
  thinning: N log-spaced points per chunk, start and end kept).
- T_ZERO = 1e-7 s initial burst added analytically on beam ON
  (`burst_dose_rates` recorded); short ON spans (< 1e-6 s at run time)
  added analytically, warning logged, recorded in `short_on_spans`;
  chunks cut exactly at beam transitions, constant `dose_rate` per chunk
  (beam at chunk midpoint), `step_dose_rates` kept for rate
  reconstruction.
- Solver failure → `SolverFailed { span, t_reached, message }` error.

### 7. Results (radiopyo `ode_result.py`)

- Block storage of chunks (O(1) append), `final_cc`.
- Rate reconstruction through the *same* `reaction_rates` function as the
  solver.
- `integrate_species`, `integrate_production_of(species)`,
  `integrate_consumption_of(species)` with radiopyo's options
  (`include_equilibrium`, `include_radiolytic`, burst production added),
  trapezoid integration, thinned-result warning.

## Semantics to preserve (radiopyo ADRs, `.work/adr/`)

The builder must read these ADRs and replicate them exactly:

- 0001 n·k convention: scalar k divided by n for `n A -> B`; polynomial,
  Arrhenius and equilibrium k used as given (warning for n > 1).
- 0002 G(H3O+) = G(e_aq).
- 0004 rate laws read clipped concentrations max(c, 0).
- 0005 smooth depletion switch c⁺/(c⁺+ε), ε = `DEPLETION_LIMIT`
  (1e-12 mol/L), for each "negative" product of k and enzymatic reactions;
  analytic Jacobian includes its derivatives.
- 0006 Michaelis–Menten k used as given.
- 0007 overshoot relaxation: a dynamic species below 0 is read as 0
  everywhere, its own row continued linearly (dc/dt += J_jj(0⁺)·c), also
  in the Jacobian.
- 0008 short ON spans added analytically.
- 0003 (as amended by 0005): analytic Jacobian is the default.

## Outputs (approach A)

Library first; a CLI binary `radkin` on top:

```
radkin run <config.toml> --t-span 1e-7 10 [--mode single|chunked|single-pulse]
           [--t-eval <file.txt>] [--rtol ..] [--atol ..] [--max-step ..]
           [--max-step-on ..] [--max-step-off ..] [--points-on N] [--points-off N]
           [--beam <label>] [--let <keV/um>] [--temperature <degC>]
           [--production O2,H2O2] [--consumption O2] --out <dir>
```

Each run writes into `<dir>/`:

| File | Mirrors | Layout / units |
|---|---|---|
| `concentrations.csv` | `res.to_pandas()` | first column `t` [s], then dynamic species (sorted order) then fixed species, µmol/L |
| `reaction_rates.csv` | `res.reconstruct_reactions()` | first column `t`, columns in radiopyo order: k-reactions (acid/base, equilibrium fwd/bwd first), radiolytic, enzymatic; labels = `as_label()`; mol/L/s |
| `species_integrals.csv` | `res.integrate_species()` | `species,value`, µmol/L·s |
| `production_<sp>.csv` / `consumption_<sp>.csv` | `integrate_production_of` / `integrate_consumption_of` | `reaction,value`, µmol/L, same sorting as radiopyo |
| `run.json` | metadata | config path, merged config, beam, LET, temperature, G-values used, burst dose rates, short ON spans, thinned flag, solver stats (steps, f/J evals), radkin version, rsoda rev |

CSV numbers written with full precision (`{:.17e}` or shortest
round-trip). The same writers are available from the library.

## Verification

1. **Golden export** — `tools/export_golden.py` (run once in the
   `radiobio` conda env) converts radiopyo's `tests/golden/data/*`
   (`runs.npz`, `derive.json`, `beams.npz`, `paper_tonneau_2025*.npz`, and
   the `.golden_results/*.json`) into JSON under `tests/golden/`. Committed.
2. **Rust tests** — unit tests per module (parser, merge, grammar, k
   values, beams/timers, G-values) ported from radiopyo's pytest suite
   where they test behaviour (not Python-specific APIs); `derive` and
   `jacobian` compared to `derive.json` at rtol 1e-12 (deterministic
   math, no solver involved); beams compared to `beams.npz` exactly;
   full runs compared to `runs.npz` / paper data at the solver tolerance
   (target rtol 1e-6 on species above a noise floor; exact thresholds
   tuned and recorded in the test, never loosened silently).
3. **Compare script** — `tools/compare_radiopyo.py <config> [beam/run
   options]`: runs radiopyo with the same inputs and the same `t_eval`,
   runs `radkin`, loads both, prints max relative error per species and per
   reaction, exit code ≠ 0 above threshold. This is the acceptance tool.

## Acceptance

Phase 1 is done when `cargo test` passes (golden data included) and
`compare_radiopyo.py` passes for: pure water, Labarbe, the paper config
(`CONFIG_PAPER_01.toml`), each with constant, single-pulse, pulsed and
macro-pulsed beams, in single and chunked modes.

## Phase-2 notes (not built now, shapes Phase 1)

- State layout must allow `n_cells × n_species` (cell-major), with the
  per-cell kernel above.
- rsoda only offers dense and banded Jacobians. A 3D grid × ~40 species is
  far beyond dense LU and does not give a useful band. Phase 2 will need
  operator splitting (rsoda per voxel for reactions + separate diffusion
  step) or a sparse/Krylov solver added to rsoda. Decided in the Phase-2
  design, not now.

## Left to the builder

- Crate/module layout, error type design (`thiserror` fine), choice of
  TOML/CSV/JSON crates (`toml`, `serde`, `serde_json`, `csv` acceptable),
  CLI crate (`clap`).
- How the scipy `solve_ivp` LSODA stepping is emulated with rsoda
  (e.g. `step_tcrit` with `tcrit = t_end` and rsoda interpolation for
  `t_eval`), as long as the acceptance tolerances are met.
- CI (GitHub Actions: `cargo fmt --check`, `clippy -D warnings`,
  `cargo test`).
