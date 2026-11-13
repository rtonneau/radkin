# radkin Phase 1 (radiopyo port) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Rust crate + CLI that reproduces radiopyo's homogeneous radiolysis results with rsoda, writing CSV files that mirror radiopyo's DataFrames.

**Architecture:** config TOML → `Config` (merged, includes resolved) → `SimEnv` (species, reactions, indices) → `SimMatrices` → per-cell kernel (`reaction_rates`, `derive`, `jacobian`) → `solve_ivp`-equivalent driver over `rsoda::Lsoda` (ITASK 5 stepping + `t_eval` interpolation) → `UnitCell` run modes → `OdeResult` → CSV/JSON writers. Each module ports one radiopyo module; the Python source (`D:\DEV\PYTHON\radiopyo\radiopyo\...`) is the reference for every behaviour not spelled out here.

**Tech Stack:** Rust 2024 stable; `rsoda` (git, rtonneau/rsoda), `toml` (feature `preserve_order`), `indexmap`, `serde`, `serde_json`, `thiserror`, `log` + `env_logger`, `clap` (derive). Python 3.11 `radiobio` env (`E:\miniconda3\envs\radiobio\python.exe`) for the tools only.

**Spec:** `docs/superpowers/specs/2026-11-13-radiopyo-port-design.md`

## Global Constraints

- Edition 2024, no C dependencies, `cargo clippy --all-targets -- -D warnings` and `cargo fmt --check` clean.
- rsoda: `rsoda = { git = "https://github.com/rtonneau/rsoda", rev = "ae76eee" }`.
- License LGPL-3.0-or-later; bundled data copied from radiopyo with attribution in `data/README.md`.
- Units: config concentrations mol/L; solver state µmol/L (dynamic species only); rates mol/L/s; Jacobian of `derive` in the µmol/L state (equals the molar Jacobian).
- Constants (exact): `T_ZERO = 1e-7`, `FIRST_STEP = 1e-15`, `MIN_BEAM_ON_TIME = 1e-6`, `_ON_TIME_RTOL = 1e-9`, `COUNT_RTOL = 1e-9`, `DEPLETION_LIMIT = 1e-12`, `DEFAULT_RTOL = 1e-10`, `DEFAULT_ATOL = 1e-12`, `DEFAULT_TRANS_TIME = 1e-4`, `THINNING_FIRST_OFFSET = 1e-6`, `GAS_CONSTANT = 8.314462618`, `ELEMENTARY_CHARGE = 1.602176634e-19`, `AVOGADRO = 6.02214076e23`, default temperature 25 °C, `0 °C = 273.15 K`, H2O `cc_init = [55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9]` (polynomial of °C), GF fit `0.97207704 * dr^-0.10972`.
- Include prefix `radiopyo/` maps to the crate `data/` dir (`radkin/` is an alias).
- Floating-point order: mass-action products multiply slots left to right, then `* k`; polynomials by numpy Horner (`c0 = c[-1]; for i in 2..=n { c0 = c[n-i] + c0*x }`); scalar k stored as `log10(k)` and evaluated as `10f64.powf(poly(1/T_K))` (k ≤ 0 → coefficient `-inf` → k = 0).
- Python float semantics where radiopyo relies on them: `//`, `%` (CPython `float_divmod`), `round(x, 6)` (via `format!("{:.6}")` round-trip).

## Review Focus

1. Reaction strings written without spaces or with `+`/`-` charges (`"H_r + H_r-> H2"`, `"Fe++ + B -> C"`, `"HO2-+H->X"`, `"2O2"`, `"A-> B - 2B + 3 C"`) must give radiopyo's labels and stoichiometry (net −1 → negative product). Test in Task 3.
2. `t_eval` containing times before `t_span[0] + T_ZERO` (dropped with a warning) or a run whose `t_span` is shorter than `T_ZERO` (error "Too Short beam pulse"). Test in Task 9.
3. A state with a dynamic species below 0 (overshoot path) in `derive` and `jacobian`; the analytic Jacobian must match finite differences there. Test in Task 7.
4. An include cycle (`a.toml` includes `b.toml` includes `a.toml`) → error naming the chain, not a stack overflow; a missing include → error with the path. Test in Task 5.
5. Solver failure (e.g. `max_step` tiny with a huge span → `MaxStepsExceeded`) returns `Err(SolverFailed{..})`, never a partial result. Test in Task 8.

---

### Task 1: Crate scaffold, numeric helpers, bundled data, CI

**Files:**
- Create: `Cargo.toml`, `LICENSE` (LGPL-3.0 text), `README.md` (stub, finished in Task 13), `.github/workflows/ci.yml`, `rustfmt.toml` (`max_width = 100`)
- Create: `src/lib.rs`, `src/error.rs`, `src/consts.rs`, `src/pyfloat.rs`, `src/poly.rs`
- Create: `data/` ← copy `radiopyo/data/{advanced_water_radiolysis.toml, reactions_biology.toml, config_Labarbe.toml}`, `radiopyo/data/G_VALUES/gvalues_from_Boscolo_2020.csv`, `examples/tonneau_2025_pmb/data/CONFIG_PAPER_01.toml`; `data/README.md` (source + LGPL attribution)

**Interfaces — Produces:**
- `pub enum RadkinError` (thiserror) with variants `Parse{input, msg}`, `Config(String)`, `Io{path, source}`, `Beam(String)`, `ShortBeamOn(String)`, `NoBeamDefined`, `SolverFailed{start, stop, t_reached, message}`, `GValue(String)`, `Value(String)`; `pub type Result<T> = std::result::Result<T, RadkinError>`.
- `consts::*` (values in Global Constraints), `consts::data_dir() -> PathBuf` (`env!("CARGO_MANIFEST_DIR")/data`, overridable by env `RADKIN_DATA_DIR`).
- `pyfloat::{floordiv(a,b) -> f64, pymod(a,b) -> f64, round6(x) -> f64}` (CPython `float_divmod`, including `x // inf == 0` for x ≥ 0).
- `poly::Polynomial { coef: Vec<f64> }` with `eval(x)`, `len()`, `from_value(v: f64, log10: bool)`, `from_coefs(Vec<f64>)` (len 1 + log10 → `log10(c0)`), `sub/add` (coefficient-wise, zero-padded), mirroring `polynomial_factory` (negative scalar → `[-inf]` when log10).

- [ ] **Step 1:** Write tests in `src/pyfloat.rs` / `src/poly.rs`: `floordiv(7.0, 2.0) == 3.0`, `floordiv(-7.0, 2.0) == -4.0`, `pymod(-7.0, 2.0) == 1.0`, `floordiv(5.0, f64::INFINITY) == 0.0`, `pymod(5.0, f64::INFINITY) == 5.0`, `floordiv(0.3, 0.1) == 2.0` (CPython gives 2.0), `round6(0.1234565) == 0.123456` / `round6(2.0000004999) == 2.0`; `Polynomial::from_coefs(vec![1.,2.,3.]).eval(2.0) == 17.0`; `from_value(1e10, true).eval(x) == 10.0`; `from_value(0.0, true)` gives k = 0 via `10^eval`.
- [ ] **Step 2:** `cargo test` → fails (missing items).
- [ ] **Step 3:** Implement. `Cargo.toml`: package `radkin` 0.1.0, `license = "LGPL-3.0-or-later"`, deps per Tech Stack, `[[bin]] name = "radkin" path = "src/bin/radkin.rs"` (stub `fn main(){}`), dev-dep `approx`. CI: windows-latest + ubuntu-latest, `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, `cargo test`.
- [ ] **Step 4:** `cargo test` → pass; `cargo clippy --all-targets -- -D warnings` clean.
- [ ] **Step 5:** Commit `feat: crate scaffold, numeric helpers, bundled data`.

### Task 2: Golden export tool

**Files:**
- Create: `tools/export_golden.py`, `tests/golden/*.json` (generated, committed), `tools/README.md`

**Interfaces — Produces** (JSON files read by Rust tests in Tasks 6–12, via `serde_json::Value`):
- `tests/golden/derive.json`: copy of radiopyo `tests/golden/data/derive.json` (keys `labarbe|pure_water` → `g_values{G_x: ge}`, `dydt{"<state>__<case>": {species: value}}`), plus per config `"states": {state: {species: y}}` recomputed with radiopyo's `_derive_states` so Rust needs no crc32 code, plus `"jacobian": {"<state>__beam_on": [[...]]}` (dense, dyn species order) and `"species_order": [labels]`, `"reaction_str": [str(r) for all k, radiolytic, enzymatic]`, `"reaction_label": [as_label()...]`.
- `tests/golden/beams.json`: from `beams.npz` — `{name: {"scalars":[t_max, dr, dr_peak], "windows":[[start,stop],...], "matrix":[[[t,dr]...]...], "intervals":[...]}}` plus the beam constructor parameters (copied from `_beam_cases`).
- `tests/golden/runs.json`: from `runs.npz` — for prefixes `labarbe_run`, `labarbe_curve`, `chunked_run`: `time`, `n_points`, `species_labels`, `species` (rows), `rate_labels` (`str(reaction)`), `rate`; `chunked_run` also `g_labels`, `g_values`; plus `"labarbe_curve_t_eval"` (the geomspace grid).
- `tests/golden/paper_points.json`: a subset (≤ 12 points: 3 dose rates × 4 pO2, both scenarios) of `paper_tonneau_2025.npz` / `paper_tonneau_2025_original.npz` with the inputs needed to rerun each point (dose rate, dose, pO2→O2 mol/L, LET, config path, `variable_GValues`, run options of `paper_runs._solve`) and the recorded outputs (final O2, H2O2, min pH or whatever keys the npz holds — exporter lists every key).
- Extra keys computed by the exporter with radiopyo: `derive.json[config]["n_reactions"]` (`{acid_base, equilibrium, k, enzymatic, radiolytic}` counts); `runs.json["thinned_times"]` (`thinned_times(0,1,5)`, `(1e-7,1e-3,10)`, `(2.0,2.0,4)`); `runs.json["labarbe_short_integrals"]` (Labarbe, `run((1e-7,1e-3), t_eval=geomspace(1e-7,1e-3,50))`: `time`, `species`, `integrate_species`, `integrate_production_of("H2O2")`, `integrate_consumption_of("O2")`, all as `{label: value}`).
- `tests/golden/parser_labels.json`: `{reaction_string: {"reactants": [[label, stoi, "+"|"-"]], "products": [...], "enzyme": label|null}}` for every reaction string in the bundled configs plus the Review-Focus edge cases.

- [ ] **Step 1:** Write `tools/export_golden.py` (argparse `--radiopyo D:\DEV\PYTHON\radiopyo --out tests/golden`; imports radiopyo from that checkout; `-I`-safe: adds the path explicitly).
- [ ] **Step 2:** Run `E:\miniconda3\envs\radiobio\python.exe -I tools/export_golden.py --radiopyo D:/DEV/PYTHON/radiopyo --out tests/golden` → five JSON files written; print their sizes.
- [ ] **Step 3:** Commit `test: golden data exported from radiopyo`.

### Task 3: Species labels and reaction-string parser

**Files:**
- Create: `src/chemistry/mod.rs`, `src/chemistry/parser.rs`
- Test: `tests/parser.rs`

**Interfaces — Produces:**
- `pub enum Operation { Add, Subtract }` (`Display`: `+`/`-`).
- `pub struct ParsedSpecies { pub label: String, pub stoi: u32, pub ope: Operation }`.
- `pub enum ReactionKind { K, Constant, Equilibrium, Enzymatic }`.
- `pub struct GenericReaction { pub kind: ReactionKind, pub reactants: IndexMap<String,(u32,Operation)>, pub products: IndexMap<String,(u32,Operation)>, pub enzyme: Option<String> }`.
- `pub fn parse_reaction(s: &str) -> Result<GenericReaction>`; `pub fn parse_species(s: &str) -> Result<String>` (label).
- `pub fn revert(g: &GenericReaction) -> Result<GenericReaction>` (error if any Subtract product).
- `pub fn reaction_key(s: &str) -> Result<ReactionKey>` where `ReactionKey = (BTreeSet<(String,u32)>, BTreeSet<(String,u32,Operation)>, Option<String>)` (for config merge).

Rules (from the Lark grammar + transformer, verified against radiopyo):
- Split on the arrow first: `<->` (Equilibrium), `--` … `>>` (Enzymatic, enzyme = plain word), `->` (K; empty left side = Constant).
- Species token: optional leading integer multiplicity (`2O2` → 2 × `O2`), then a formula or a plain element. Formula = one or more `Capital[lower]?[digits]?`, then optional `_r`, then charges `_plus`/`_minus`/`+`/`-` repeated; label = atoms joined (stoi printed when > 1) + `_r` + `_plus`×n / `_minus`×n; `e_aq` is its own label. Anything else matching `\w+` is an element whose label is the text.
- Attached `+`/`-` run after a token: if the next non-space char starts a species, the last char of the run is the separator, the rest are charges; otherwise the whole run is charge.
- Reactants: duplicates summed (`H_r + H_r` → `H_r`×2). Products: signed sums; net < 0 → `(abs, Subtract)`; net 0 removed. Enzymatic: exactly one reactant.

- [ ] **Step 1:** Write `tests/parser.rs`: load `tests/golden/parser_labels.json`, and for each entry assert `parse_reaction` gives the same reactants/products (label, stoi, op) in the same order and the same enzyme; plus explicit asserts: `"A-> B - 2B + 3 C"` → products `[("B",1,Subtract),("C",3,Add)]`; `"2H2O2 -- catalase >> O2 + 2 H2O"` → enzyme `catalase`, reactants `[("H2O2",2,Add)]`; `"-> A + B"` → kind Constant; `"A + B -> "` parses with `products.is_empty()` (radiopyo accepts it); `"A + -> B"` is `Err(Parse)` whose message contains the string; `"A -- E >> B + C - D"` ok; `parse_species("OH_r")=="OH_r"`, `parse_species("Fe++")=="Fe_plus_plus"`.
- [ ] **Step 2:** `cargo test --test parser` → fails.
- [ ] **Step 3:** Implement a hand-written tokenizer + parser in `src/chemistry/parser.rs`.
- [ ] **Step 4:** `cargo test --test parser` → pass.
- [ ] **Step 5:** Commit `feat: reaction-string parser with radiopyo labels`.

### Task 4: Timers and beams

**Files:**
- Create: `src/physics/mod.rs`, `src/physics/timer.rs`, `src/physics/beam.rs`, `src/physics/beam_collection.rs`
- Test: `tests/beams.rs`

**Interfaces — Produces:**
- `pub enum TimerState { Off, On }`; `pub struct Timer { period, on_time, start_time }` with `state_at`, `time_till_next_on/off`, `time_till_next_period`, `period_index` (via `pyfloat::floordiv`), `duty_cycle`, `is_continuous`; `pub struct MacroTimer { period, on_time, start_time, sub: Timer }` with `state_at`, `micro_pulses_per_train`, `duty_cycle`, `time_till_next_on/off`. `split_periods(duration, period) -> (i64, f64)`.
- `pub enum BeamKind { Constant, SinglePulse, Pulsed, MacroPulsed { macro_timer: MacroTimer } }`; `pub struct Beam { pub kind, pub dose_rate_peak, pub timer: Timer, pub max_dose, pub label, pub time_at_max_dose, pub dose_rate }`.
- Constructors: `Beam::constant(dose_rate, max_dose: Option<f64>, start_time: Option<f64>) -> Result<Beam>`, `Beam::single_pulse(dose_rate, max_dose, start_time)`, `Beam::single_pulse_from_on_time(dose_rate, on_time, start_time)`, `Beam::pulsed(dose_rate_peak, period, on_time, max_dose, start_time)`, `Beam::pulsed_from_dose_rate_peak(peak, period, on_time, n_pulse: Option<f64>, max_dose: Option<f64>, start_time)` (isclose 1e-9 consistency check), `Beam::pulsed_from_dose_per_pulse(..)`, `Beam::macro_pulsed(peak, period, on_time, micro_period, micro_on_time, max_dose, micro_start_time, start_time)`; all run `post_init` (dose_rate_peak ≤ 0 error, `time_at_max_dose`, average `dose_rate`, ShortBeamOn check).
- `Beam::from_params(beam_type: &str, params: &IndexMap<String, f64>) -> Result<Beam>` (radiopyo `from_dict` of each type, all parametrisations), `Beam::as_params(&self) -> IndexMap<String,f64>` (radiopyo `as_dict`), `Beam::beam_type() -> &str` (`"constant"`, `"singlepulse"`, `"pulsed"`, `"macro_pulsed"`, lower-case as in `as_dict` merge).
- `Beam::at(t) -> f64` (dose rate), `Beam::tot_pulses()`, `Beam::dose_per_pulse()`, `Beam::make_matrix(stop, start) -> Vec<[f64;2]>`, `Beam::make_time_intervals(stop, start) -> Vec<[f64;2]>` (second column 0/1).
- `pub struct BeamCollection` (IndexMap label → (Beam, Option<LET>), `current`): `add_beam(beam, key, let_, use_)`, `use_beam`, `current()`, `current_let() -> Result<f64>`, `update_let`, `remove`, `from_config(&IndexMap<String, Vec<IndexMap<String, toml::Value>>>) -> Result<Self>` (pops `LET`, `label` default `"default"`, strings other than label rejected), `as_config()`, `merge(&self, other) -> Result<Self>` (radiopyo `BeamCollection.merge`).

- [ ] **Step 1:** Write `tests/beams.rs`: for every case in `tests/golden/beams.json`, rebuild the beam from its parameters, assert `[time_at_max_dose, dose_rate, dose_rate_peak]` and every `make_matrix` / `make_time_intervals` window equal golden with rel 1e-12; plus `Beam::pulsed(1e6, 0.01, 5e-7, None, None)` → `Err(ShortBeamOn)`; `from_params("pulsed", {frequency: 90, on_time: 2e-6, dose_rate_peak: 5.7e5, max_dose: 10})` period = 1/90; `merge` of a `default` constant with a labelled pulsed keeps both.
- [ ] **Step 2:** `cargo test --test beams` → fails.
- [ ] **Step 3:** Implement (port `timer.py`, `beam.py`, `beam_collection.py`; transitions as an iterator struct or a loop that builds rows directly).
- [ ] **Step 4:** `cargo test --test beams` → pass.
- [ ] **Step 5:** Commit `feat: timers, beams, beam collection`.

### Task 5: Config loading, includes, merge

**Files:**
- Create: `src/config/mod.rs`, `src/config/merge.rs`, `src/config/include.rs`
- Test: `tests/config.rs`, fixtures `tests/fixtures/config/{base.toml, override.toml, cycle_a.toml, cycle_b.toml, sub/nested.toml}`

**Interfaces — Produces:**
- `pub struct Config { pub includes: IndexMap<String, Vec<String>>, pub parameters: toml::Table, pub fixed: IndexMap<String, toml::Value>, pub initial: IndexMap<String, toml::Value>, pub radiolytic: IndexMap<String, f64>, pub reactions: IndexMap<Section, Vec<toml::Table>>, pub beam: IndexMap<String, Vec<toml::Table>> }`, `pub enum Section { AcidBase, Equilibrium, KReaction, Constant, Arrhenius, Enzymatic }` (TOML names `acid_base`, `equilibrium`, `k_reaction`, `constant`, `arrhenius`, `enzymatic`).
- `Config::from_toml_file(path) -> Result<Config>` (parse + `anchor_includes`), `Config::load(path) -> Result<Config>` (= `resolve_includes(from_toml_file(path))`), `combine(root, other) -> Result<Config>` (radiopyo `combine_config`, including the beam merge through `BeamCollection` and the `default`-beam drop), `resolve_includes(cfg, chain) -> Result<Config>`.
- `Config::parameters` accessors: `ph() -> f64` (7), `force_cst_ph() -> bool`, `temperature() -> f64` (25), `variable_gvalues() -> bool`, `species_order() -> Option<(Vec<String>, Vec<String>)>`.

- [ ] **Step 1:** Write `tests/config.rs`: `Config::load(data/config_Labarbe.toml)` has 1 constant beam (`dose_rate=1e7`, LET 0.201) and the k-reaction count of both includes (radiopyo: 64+13 strings minus overrides — assert against a count exported in `tests/golden/derive.json["labarbe"]["n_reactions"]` added in Task 2's exporter); override fixture: `"B + A -> C"` in override replaces `"A + B -> C"` k_value in place; equilibrium override giving `k_forward` + `k_equilibrium` drops old `k_backward`; acid_base override with only `k_forward` keeps `pKa`; nested include resolved relative to including file; cycle fixture → `Err` message containing both file names; missing include → `Err` containing the path; `radiopyo/` prefix resolves to `data/`.
- [ ] **Step 2:** `cargo test --test config` → fails.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `cargo test --test config` → pass.
- [ ] **Step 5:** Commit `feat: TOML config, includes and merge rules`.

### Task 6: Chemistry objects, SimEnv, G-values

**Files:**
- Create: `src/chemistry/reactions.rs`, `src/chemistry/species.rs`, `src/chemistry/env.rs`, `src/gvalues.rs`
- Test: `tests/env.rs`, unit tests in `src/gvalues.rs`

**Interfaces — Produces:**
- `pub struct SimSpecies { pub label: String, pub index: usize, pub cc_init: Option<Polynomial>, pub fixed: bool }` with `initial_cc(t_celsius) -> f64` (clamped ≥ 0).
- `pub struct ListSpecies` (IndexMap dyn, IndexMap cst; `get_or_create_dyn`, `add_dyn`, `add_cst`, `sort(starter)`, `set_order(dyn, cst)`, `n_dyn()`, `n_cst()`, `sorted_dyn()`, `sorted_cst()`, `is_cst(label)`, `index_of(label)`).
- `pub struct ReactionSpecies { pub label: String, pub stoi: u32, pub ope: Operation }`.
- `pub enum KValue { Poly(Polynomial) /* log10 k of 1/T_K */, Arrhenius { a: f64, ea_kj: f64 } }`; `pub struct KReaction { pub reactants: Vec<ReactionSpecies>, pub products: Vec<ReactionSpecies>, pub k: Option<KValue>, pub index: usize, pub kind: KKind /* K, Constant, Arrhenius, Forward, Backward */ }` with `k_value(t_celsius) -> Result<f64>`, `as_label()` (`KR_i` / `CR_i`), `Display` = radiopyo `__str__` (`as_string` of each side joined by `->` or `<->`; Constant: `"-> " + products`).
- `pub struct BiReaction { pub forward: KReaction, pub backward: KReaction, pub pka: Option<Polynomial> /* acid/base */, pub acid: Option<String>, pub base: Option<String> }` (acid/base couple detection per `AcidBase.get_couple`; equilibrium constants per `EReaction.__init__`).
- `pub struct MichaelisMenten { substrate stoi, reactants, products, enzyme: String, k_value, k_micha, index }` (`EZ_i`, `Display` `"A -- E >> B"`).
- `pub struct RadiolyticReaction { pub species: String, pub kr: f64, pub index: usize }` (`ge_to_kr(ge) = ge / ELEMENTARY_CHARGE / 100 / AVOGADRO`, `kr_to_ge`, `as_label() = "G_<sp>"`, `Display = "RadiolyticReaction('<sp>')"`).
- `pub struct ReactionsList { pub ab: Vec<BiReaction>, pub eq: Vec<BiReaction>, pub k: Vec<KReaction> /* k_reaction, then arrhenius, then constant */, pub enzymatic: Vec<MichaelisMenten>, pub radiolytic: Vec<RadiolyticReaction> }`; `iter_all_k()` (ab fwd,bwd…, eq fwd,bwd…, k…), global indices per radiopyo order, `involving(species, role, kinds)` queries needed by Task 10.
- `pub struct SimEnv { pub reactions, pub species, pub temperature: f64, pub ph: f64, pub force_cst_ph: bool, pub variable_gvalues: bool }`; `SimEnv::from_config(&Config) -> Result<SimEnv>` (radiopyo `_chemistry_builder`: fixed species first, H2O fixed, H3O_plus/OH_minus dyn or fixed, reactions in section order, n·k rule, initial concentrations, `trash` → fixed, ordering by acid/base starter or `species_order`), `initial_values() -> Vec<f64>` (µmol/L, dyn), `set_temperature(c)`.
- `pub struct GValuesBoscolo` with `load(data_dir) -> Result<Self>`, `contains(species)`, `value_at(species, let_, dose_rate: Option<f64>) -> Result<f64>` (H3O_plus → e_aq; LET outside [x0, xN] → `Err(GValue)`); quadratic interpolant = scipy `make_interp_spline(x, y, k=2)`: knots `t = [x0,x0,x0, (x[i]+x[i+1])/2 for i in 1..n-2, xN,xN,xN]`, solve the n×n collocation system (dense LU with partial pivoting), evaluate by de Boor.

- [ ] **Step 1:** Write `tests/env.rs`: for `labarbe` and `pure_water` (cell built as in radiopyo `build_cell`: config's beam, else `Beam::constant(1e7, Some(10.0), None)` with LET 0.201), assert species order == golden `species_order`, `Display` of all reactions == golden `reaction_str`, labels == golden `reaction_label`; `KReaction` for `"2 OH_r -> H2O2"` with k_value 1.1e10 gives `k_value(25) == 5.5e9` (rel 1e-14); Arrhenius `A=1e10, Ea=20` → `3.13e6` within 1%; unit tests in `gvalues.rs`: `value_at("e_aq", 0.13, None) == 2.7` (rel 1e-12), `value_at("H3O_plus", x, None) == value_at("e_aq", x, None)`, `value_at("OH_r", 200.0, None)` is `Err`; golden check: after `set LET 0.201`, the G-values `g_values` in `derive.json` match `kr_to_ge(kr)` rel 1e-12.
- [ ] **Step 2:** `cargo test --test env` → fails.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `cargo test --test env` and `cargo test --lib gvalues` → pass.
- [ ] **Step 5:** Commit `feat: chemistry objects, SimEnv builder, Boscolo G-values`.

### Task 7: Matrices and per-cell kernel (rates, derive, Jacobian)

**Files:**
- Create: `src/matrices.rs`, `src/kernel.rs`
- Test: `tests/kernel.rs`

**Interfaces — Consumes:** `SimEnv` (Task 6).
**Produces:**
- `pub struct KMatrices { k_values: Vec<f64>, reactant_index: Vec<Vec<usize>> /* padded with n_species */, stoichiometry: Vec<f64> /* n_species × n_kr row-major */, negative_products: Vec<bool>, negative_index: Vec<Vec<usize>> }`, `EzMatrices { k_values, k_michaelis, substrate_index, enzyme_index, stoichiometry, negative_products, negative_index }`, `GVector { yields, stoichiometry }`, `StateLayout { cst_sp: Vec<f64> }`, `OvershootTables` (per dyn species `Option<KOvershootTable>`, `Option<EzOvershootTable>`, fields as radiopyo `sim_types.py`), `pub struct SimMatrices { k, radiolytic, enzymatic, state, overshoot, n_species, n_dyn }`.
- `SimEnv::build_matrices(&self) -> SimMatrices` (radiopyo `build_all_matrices`).
- `pub struct Workspace` (preallocated buffers: `y_full`, rates, switch, tmp) — `Workspace::new(&SimMatrices)`.
- `pub fn reaction_rates(y_full: &[f64], dose_rate: f64, mat: &SimMatrices, out: &mut ReactionRates)` with `pub struct ReactionRates { pub k: Vec<f64>, pub radiolytic: Vec<f64>, pub enzymatic: Vec<f64> }` (`stacked()` → k, radiolytic, enzymatic).
- `pub fn derive(y: &[f64], dose_rate: f64, mat: &SimMatrices, ws: &mut Workspace, dydt: &mut [f64])` (µmol/L in and out, overshoot relaxation).
- `pub fn jacobian(y: &[f64], mat: &SimMatrices, ws: &mut Workspace, pd: &mut [f64])` — dense `n_dyn × n_dyn`, column-major (`pd[i + j*n_dyn] = ∂f_i/∂y_j`, rsoda layout). Dose rate does not enter the Jacobian.
- These take one cell's slices and allocate nothing (Phase-2 seam).

- [ ] **Step 1:** Write `tests/kernel.rs`: for each config × state (`initial`, `synthetic`, `zeros`) × case (`beam_on` dose rate = beam.at(5e-7), `beam_off` = beam.at(1e-3), `no_beam` = 0), `derive` equals golden `dydt` with `rtol 1e-12, atol 1e-12·max|golden|`; `jacobian` at `synthetic__beam_on` equals golden `jacobian` (same tolerance); FD check: at a state with `OH_r = -1e-3` µmol/L (overshoot path) and at `synthetic`, central differences of `derive` (h = 1e-7·max(|y_j|, 1e-6)) match `jacobian` within rtol 1e-5 / atol 1e-6·max|J|.
- [ ] **Step 2:** `cargo test --test kernel` → fails.
- [ ] **Step 3:** Implement (straight loops over the padded index tables; products left to right; matmul as Σ over reactions in index order).
- [ ] **Step 4:** `cargo test --test kernel` → pass.
- [ ] **Step 5:** Commit `feat: reaction matrices and per-cell derive/jacobian kernel`.

### Task 8: `solve_ivp`-equivalent LSODA driver

**Files:**
- Create: `src/solver.rs`
- Test: unit tests in `src/solver.rs`

**Interfaces — Consumes:** `kernel::{derive, jacobian}`, rsoda.
**Produces:**
- `pub enum DoseRateSource<'a> { None, Constant(f64), Beam(&'a Beam) }` with `at(t) -> f64`.
- `pub struct SolveOptions { pub rtol: f64, pub atol: f64, pub first_step: f64, pub max_step: Option<f64>, pub use_jac: bool, pub t_eval: Option<Vec<f64>> }`.
- `pub struct Solution { pub t: Vec<f64>, pub y: Vec<Vec<f64>> /* per time point, n_dyn */, pub nfev: usize, pub njev: usize, pub nst: usize }`.
- `pub fn solve(mat: &SimMatrices, source: DoseRateSource, t_span: (f64, f64), y0: &[f64], opts: &SolveOptions) -> Result<Solution>`.

Algorithm (scipy `solve_ivp` + `LSODA` wrapper semantics):
- rsoda options: `Tolerance::scalar(rtol, atol)`, `h0 = first_step`, `hmax = max_step` (None → default), `JacobianKind::DenseUser` if `use_jac` else `DenseFd`, `mxstep` default.
- Loop `t = solver.step_tcrit(t_end, t_end)` until `t >= t_end`; keep `t_old`. Without `t_eval`: store `t0, y0` then every `(t, y)`. With `t_eval` (must lie in `[t0, t_end]`, else `Err(Value)`): after each step, the points `t_eval[i..j]` with `t_eval <= t` (searchsorted right) are interpolated by `solver.integrate(te)` (rsoda interpolates inside the last step without stepping; DINTDY = scipy `LsodaDenseOutput`), then continue stepping. `nfev = solver.nfe()`, `njev = solver.nje()`.
- Any rsoda error → `Err(SolverFailed { start, stop, t_reached: solver.t(), message: e.to_string() })`.

- [ ] **Step 1:** Unit tests: exponential decay system built from a tiny config string (`"A -> B"`, k=1, initial A=1e-3 mol/L) solved on (0, 1) with defaults: final A = 1000·e⁻¹ µmol/L rel 1e-8; same with `t_eval = [0, 0.5, 1]` gives exactly 3 points at those times, values rel 1e-8; `t_eval` with 2.0 → `Err`; `rtol = 0, atol = 0` with a species at 0 → rsoda `ErrorWeightZero` → `Err(SolverFailed)` carrying the span and the rsoda message (no partial `Solution`).
- [ ] **Step 2:** `cargo test --lib solver` → fails.
- [ ] **Step 3:** Implement `OdeSystem` adapter struct `ChemSystem<'a> { mat, source, ws }` (rhs → `derive(y, source.at(t), ..)`, jac → `jacobian`).
- [ ] **Step 4:** `cargo test --lib solver` → pass.
- [ ] **Step 5:** Commit `feat: solve_ivp-equivalent LSODA driver on rsoda`.

### Task 9: UnitCell and run modes

**Files:**
- Create: `src/unit_cell.rs`
- Test: `tests/unit_cell.rs`

**Interfaces — Consumes:** Tasks 4–8. **Produces:**
- `pub struct RunOptions { pub atol: Option<f64>, pub rtol: Option<f64>, pub max_step: Option<f64>, pub first_step: Option<f64>, pub use_beam: bool /* true */, pub use_jac: bool /* true */, pub t_eval: Option<Vec<f64>> }` (`Default`).
- `pub struct ChunkOptions { pub run: RunOptions, pub max_step_size_on, pub max_step_size_off, pub max_step_size_trans, pub max_step_size_final: Option<f64>, pub trans_time: Option<f64>, pub points_per_chunk_on, pub points_per_chunk_off: Option<usize>, pub y0: Option<Vec<f64>> }`.
- `pub struct SimChunk { pub start, pub stop, pub max_step_size: Option<f64>, pub beam_on: bool, pub dose_rate: Option<f64>, pub n_points: Option<usize> }`.
- `pub struct UnitCell { pub env: SimEnv, pub mat: SimMatrices, pub beams: BeamCollection, pub gvalues: GValuesBoscolo }`: `from_config(&Config)`, `from_toml(path)`, `set_beam(Beam, Option<f64> LET)`, `use_beam(label, Option<LET>)`, `change_let(let_)` (rewrites G of e_aq, OH_r, H_r, H2, H2O2, H3O_plus when dynamic; GF scaling with `dose_rate_peak` when `variable_GValues`), `set_temperature(c)`, `run(t_span, &RunOptions, y0: Option<&[f64]>) -> Result<OdeResult>`, `run_till_max_dose(time_start=1e-9, ..)`, `prepare_chunked_run(t_span, &ChunkOptions) -> Result<Vec<SimChunk>>`, `prepare_single_pulse(t_span, &ChunkOptions) -> Result<Vec<SimChunk>>`, `run_chunks(t_span, chunks, &ChunkOptions) -> Result<OdeResult>`.
- Free fns: `thinned_times(start, stop, n) -> Vec<f64>`, `widen_empty_span`, `is_short_on_span(span, min_on)`, `chunk_dose_rate(beam, start, stop)`, `check_points_per_chunk`.
- `_raw_run` logic exactly as radiopyo (burst at `t_span.0 + T_ZERO`, `last_burst_dose_rate`, short ON span step, `_drop_skipped_t_eval`, thinning replaces `t_eval`, `_solver_kwargs` first-step rules).

- [ ] **Step 1:** Write `tests/unit_cell.rs`: `thinned_times(0.0, 1.0, 5)` == radiopyo values (export 3 cases in Task 2's exporter under `runs.json["thinned_times"]`); `prepare_chunked_run((0, 0.05))` for the golden chunked beam gives the boundaries of golden `chunked_run` edges; `run((1e-7, 1e-7 + 5e-8))` with a constant beam → `Err` "Too Short beam pulse"; `t_eval = [1e-7, 1.5e-7, 1e-3]` on `run((1e-7, 1e-3))` drops `1e-7` (2 points kept); a pulsed beam whose last pulse is cut to 0.5 µs by `max_dose` records one `short_on_span`; a cell with no beam → `Err(NoBeamDefined)`.
- [ ] **Step 2:** `cargo test --test unit_cell` → fails.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `cargo test --test unit_cell` → pass.
- [ ] **Step 5:** Commit `feat: UnitCell run, chunked and single-pulse modes`.

### Task 10: OdeResult and analysis

**Files:**
- Create: `src/result.rs`
- Test: `tests/result.rs`

**Interfaces — Produces:**
- `pub struct ShortOnSpan { pub start, pub stop, pub dose_rate }`.
- `pub struct OdeResult { pub env: SimEnv (clone), pub mat: SimMatrices (clone), pub beam: Beam, pub let_: Option<f64>, pub t_zero_pulse: f64, t_blocks, y_blocks, pub burst_dose_rates: Option<Vec<f64>>, pub step_dose_rates: Vec<Option<f64>>, pub short_on_spans: Vec<ShortOnSpan>, pub thinned: bool, pub n_points: Vec<usize>, pub nfev: Vec<usize>, pub njev: Vec<usize>, pub exec_time: Vec<f64> }`.
- `append(sol: Solution, exec_time, bursts: &[f64], dose_rate: Option<f64>, short: Option<ShortOnSpan>, thinned: bool)`, `t() -> &[f64]`, `y() -> &[Vec<f64>]`, `final_cc() -> Vec<f64>`.
- `concentrations() -> (Vec<String>, Vec<Vec<f64>>)` (= `to_pandas`: dyn sorted, then cst sorted, µmol/L).
- `reaction_rates() -> (Vec<String> /* as_label */, Vec<Vec<f64>>)` (= `reconstruct_reactions`, per-block `ConstantDoseRate` when step dose rates are set, otherwise the beam at each t; columns: all k (ab, eq, k), radiolytic, enzymatic).
- `integrate_species(start=0, stop=inf) -> Vec<(String, f64)>`, `integrate_production_of(sp, start, stop, include_equilibrium=true, include_radiolytic=true) -> Vec<(String, f64)>` (sorted descending; burst + short-span radiolytic additions; `_burst_dose_rates` fallback), `integrate_consumption_of(sp, start, stop, include_equilibrium=true) -> Vec<(String, f64)>`; trapezoid = numpy `trapezoid` (Σ (x[i+1]-x[i])·(y[i]+y[i+1])/2); time slicing inclusive both ends; thinned warning once.

- [ ] **Step 1:** Write `tests/result.rs`: hand-built `OdeResult` on the decay system: `integrate_species` of A over [0,1] ≈ 1000·(1−e⁻¹) (trapezoid tolerance of the stored points), production of `B` == consumption of `A` (rel 1e-12); for `labarbe` short run (1e-7, 1e-3) with radiopyo equivalence numbers exported in Task 2 (`runs.json["labarbe_short_integrals"]`: `integrate_production_of("H2O2")`, `integrate_consumption_of("O2")`, `integrate_species()` on `t_eval = geomspace(1e-7, 1e-3, 50)`) rel 1e-6.
- [ ] **Step 2:** `cargo test --test result` → fails.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `cargo test --test result` → pass.
- [ ] **Step 5:** Commit `feat: OdeResult, rate reconstruction and integrals`.

### Task 11: Output writers and CLI

**Files:**
- Create: `src/output.rs`, `src/bin/radkin.rs`
- Test: `tests/cli.rs`

**Interfaces — Produces:**
- `pub fn write_run_dir(dir: &Path, res: &OdeResult, meta: &RunMeta, production: &[String], consumption: &[String]) -> Result<()>` writing `concentrations.csv` (header `t,<labels>`), `reaction_rates.csv` (header `t,<as_label…>`), `species_integrals.csv` (`species,value`), `production_<sp>.csv` / `consumption_<sp>.csv` (`reaction,value`), `run.json`; numbers formatted with `{:e}` shortest round-trip (Rust `{:e}` on f64 is round-trip).
- `pub struct RunMeta { config, mode, t_span, beam_label, beam_params, let_, temperature, g_values: IndexMap<String,f64>, burst_dose_rates, short_on_spans, thinned, n_points, nfev, njev, radkin_version, rsoda_rev, reaction_strings: IndexMap<label, str(reaction)> }` (serde).
- CLI `radkin run <config> --t-span <a> <b> [--mode run|chunked|single-pulse] [--t-eval FILE] [--rtol --atol --max-step --first-step] [--max-step-on --max-step-off --max-step-trans --max-step-final --trans-time] [--points-on N --points-off N] [--beam LABEL] [--let X] [--temperature C] [--no-jac] [--production A,B] [--consumption A,B] --out DIR`; `--t-eval FILE` = one time per line; `RUST_LOG` controls logging (default `warn`).

- [ ] **Step 1:** `tests/cli.rs` (uses `env!("CARGO_BIN_EXE_radkin")`): run `data/advanced_water_radiolysis.toml` with no beam in config → exit ≠ 0 with "No beam"; run `data/config_Labarbe.toml --t-span 1e-7 1e-3 --production H2O2 --consumption O2 --out <tmp>` → the 6 files exist, `concentrations.csv` header starts with `t,` and has dyn+cst species count + 1 columns, `reaction_rates.csv` header ends with the enzymatic labels, `run.json` parses and has `g_values`.
- [ ] **Step 2:** `cargo test --test cli` → fails.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `cargo test --test cli` → pass.
- [ ] **Step 5:** Commit `feat: CSV/JSON outputs and radkin CLI`.

### Task 12: Golden run reproduction tests

**Files:**
- Create: `tests/golden_runs.rs`

- [ ] **Step 1:** Write tests (radiopyo `test_golden.py` semantics: species `rtol 1e-6, atol 1e-10`; rates `rtol 1e-6, atol 1e-7·max|rate|`; times at edges `rtol 1e-12`; `n_points` within 15 %): `labarbe_run` (`run((1e-7, 10))`, compare first/last rows), `labarbe_curve` (`t_eval` = golden grid), `chunked_run` (pure water, `Beam::pulsed_from_dose_per_pulse(2.0, 0.01, 2e-6, Some(3.0), None, None)`, LET 2.5, `max_step_size_on=1e-7, off=1e-3`, chunk edge rows, G-values rel 1e-12); paper points from `paper_points.json` (each rerun with its recorded inputs, outputs rel 1e-4 — the paper golden tolerance used by radiopyo's `test_paper_tonneau_2025.py`, read from that file when exporting). Rate rows are matched by `str(reaction)` via `RunMeta.reaction_strings`.
- [ ] **Step 2:** `cargo test --release --test golden_runs` → investigate any failure with superpowers:systematic-debugging; never loosen a tolerance without recording why in the test.
- [ ] **Step 3:** Commit `test: golden run reproduction against radiopyo`.

### Task 13: compare_radiopyo.py, acceptance matrix, README

**Files:**
- Create: `tools/compare_radiopyo.py`, `tools/acceptance.ps1`
- Modify: `README.md`

- [ ] **Step 1:** `compare_radiopyo.py <config.toml> --mode run|chunked|single-pulse --t-span A B [--beam-type T --beam-param k=v ...] [--let X] [--t-eval-n N] [radkin CLI options] [--radkin target/release/radkin] [--rtol 1e-6]`: builds the radiopyo `UnitCell` with the same beam/options, makes a shared `t_eval` (geomspace N points inside each chunk for chunked modes is not possible through `prepare_chunked_run`, so for chunked modes it compares on radiopyo's chunk edge rows and on `points_per_chunk_*` thinned grids passed to both codes), runs radkin via subprocess into a temp dir, loads CSVs, prints max rel error per species / reaction and integrals, exits 1 above `--rtol` (with an absolute floor `1e-10` µmol/L for species and `1e-7·max` for rates).
- [ ] **Step 2:** `tools/acceptance.ps1` runs compare for {pure water, Labarbe, CONFIG_PAPER_01} × {constant, single-pulse, pulsed, macro-pulsed} × {run, chunked} (single-pulse beam also with `single-pulse` mode); summary table; exit ≠ 0 on any failure.
- [ ] **Step 3:** Run `cargo build --release` then `pwsh tools/acceptance.ps1` → all PASS (fix discrepancies via systematic debugging).
- [ ] **Step 4:** README: what radkin is, build, CLI usage, output files ↔ radiopyo DataFrames table, compare workflow, Phase-2 note, license.
- [ ] **Step 5:** Commit `feat: radiopyo comparison tool and acceptance matrix`, push.
