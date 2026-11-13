# Session: rsoda-dense-output

**Date:** 2026-11-13T07:10:33.735Z
**Status:** Grill phase complete

## Problem Statement

radkin fills `t_eval` points by calling `Lsoda::integrate(te)` inside the last step, because rsoda exposes no dense-output API. That relies on `integrate` not stepping, overwrites the solver's `t()`/`y()`, and uses DINTDY, while scipy `solve_ivp(method="LSODA")` (radiopyo's solver) uses its own `LsodaDenseOutput` evaluator. Result: a less explicit contract and interpolation differences vs radiopyo beyond rounding.

## Context & Constraints

- **Current behavior:** `src/solver.rs` steps with `step_tcrit(t_end, t_end)`; for each `t_eval` point inside the step it calls `solver.integrate(tout)` (DINTDY over NQCUR columns, Horner form, scale HCUR). `rsoda::interpolate::intdy` is `pub` but `Lsoda` exposes no Nordsieck history.
- **Pain point:** scipy's `LsodaDenseOutput` differs: order NQU (last successful order, `iwork[13]`), `h = HCUR` (`rwork[11]`), last column rescaled by `(HCUR/HU)^NQU` when NQCUR < NQU, power-series evaluation `yh · ((t - t_step)/h)^p`. With a decreasing order this is a truncation-level difference, not only rounding.
- **Dependencies:** rsoda pinned by git `rev = "ae76eee"` in `Cargo.toml`; must not switch to a `path` dependency in committed code. Golden tolerances (1e-12 kernel/G-values, rtol 1e-6 runs) must not be loosened. Existing `nqu()`, `nqcur()`, `hu()`, `hcur()`, `tcur()` accessors on `Lsoda`.
- **Tech stack:** Rust 2024, rsoda (LSODA port), scipy `scipy/integrate/_ivp/lsoda.py` as reference (`LsodaDenseOutput`, `_dense_output_impl`, `solve_ivp` t_eval loop).

## Success Metrics

- rsoda: `Lsoda::dky(t, k, &mut out)` (faithful DINTDY, errors before the first step and for illegal t/k) and a read-only Nordsieck view (`yh`, `nyh`, `n`, step data) are public, unit-tested, merged to rsoda main via PR.
- radkin: `rev` bumped to the rsoda merge commit; `solve()` fills `t_eval` with an exact port of scipy's `LsodaDenseOutput` and no longer calls `integrate(te)`; a golden test checks the evaluator against scipy's `LsodaDenseOutput` on fixed `yh` arrays (including the order-decrease rescale case) at 1e-12.
- `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` clean in both repos; the 27-case acceptance matrix passes; `compare_radiopyo` on a `t_eval` case reported before/after.

## Architecture & Approach

1. rsoda (branch + PR): add `Lsoda::dky(t, k, out) -> Result<(), LsodaError>` wrapping `intdy` with the current `StepData` (mapped like `interpolate`, `InputError::Interpolation`), and `Lsoda::nordsieck()` returning a borrowed view of the history array and its shape/step data. No stepping, no change to `t()`/`y()`. No version bump (git dependency).
2. radkin: bump `rev`. Add a scipy dense-output port (in `src/solver.rs` or a small module): build `yh` columns `0..=NQU` from the Nordsieck view, rescale the last column by `(HCUR/HU)^NQU` when NQCUR < NQU, evaluate `sum_j yh[:, j] * x^j` with `x = (t - t_step)/HCUR`, `x^j` via `powf` (numpy `**`). The reference time `t_step` is the value returned by `step_tcrit` (equals `tcrit` exactly at the span end, as scipy's `self.t`).
3. `solve()` uses the evaluator for every `t_eval` point in `(t_prev, t_step]`; the `integrate(te)` path and its doc comment go.
4. Golden: `tools/export_golden.py` exports scipy `LsodaDenseOutput` outputs for synthetic `yh` cases; a Rust test compares at 1e-12.

## Assumptions & Trade-offs

- rsoda's history after a step matches ODEPACK's (rescaled to HCUR for the next step), since rsoda is a faithful port; scipy's reconstruction therefore applies unchanged.
- `np.dot` goes through BLAS with an unspecified summation order; bitwise equality with scipy is not claimed, only rounding-level agreement.
- The scipy-specific evaluator lives in radkin, not rsoda; rsoda stays an LSODA port and exposes only DINTDY plus raw history.
- Not doing: sparse Jacobian, per-step overhead work (separate ideas).

## Open Questions

None.

## Notes

- Promotes idea `.work/idea/rsoda-dense-output.md`.
- Design approved; user chose `/gps auto` to run plan, ship and finish without further stops.
