# Implementation Plan

**Session:** rsoda-dense-output
**Date:** 2026-11-13T07:13:19.046Z
**Estimated effort:** 0.5 day

## Strategy

Land the rsoda API first (its own branch and PR in rtonneau/rsoda, merged to rsoda main), because radkin can only pin a pushed revision. In parallel-safe order in radkin: first the pure scipy dense-output evaluator with its golden test (no rsoda dependency, works on slices), then the rev bump that wires it into `solve()` through the new Nordsieck view and removes the `integrate(te)` path, verified by the acceptance matrix.

## Tickets Overview

- **Ticket 1:** rsoda: public `Lsoda::dky` (DINTDY) and read-only `Lsoda::nordsieck()` view, tested, merged via PR.
- **Ticket 2:** radkin: exact port of scipy `LsodaDenseOutput` as a pure function, golden-tested against scipy at 1e-12.
- **Ticket 3:** radkin: bump rsoda rev, use the evaluator for `t_eval` in `solve()`, drop `integrate(te)`; acceptance matrix and before/after compare.

## Sequencing Rationale

Ticket 3 needs ticket 1 merged (new rev) and ticket 2's evaluator. Ticket 2 is independent of rsoda and keeps radkin checks green on its own.

## Risks & Mitigation

- **Risk:** rsoda's post-step history differs from ODEPACK's layout that scipy assumes (columns rescaled to HCUR, NQU vs NQCUR) → **Mitigation:** unit test in rsoda that `dky(tn, 0)` equals `y()` at a mesh point and that history columns match `intdy`; radkin acceptance matrix at rtol 1e-6.
- **Risk:** BLAS summation order in `np.dot` makes bitwise comparison impossible → **Mitigation:** golden tolerance 1e-12 relative, as the other goldens; no bitwise claim.
- **Risk:** PR merge on rsoda blocked (CI) → **Mitigation:** fix in the rsoda branch; radkin ticket 3 waits on the merge commit.

## Assumptions

- `${DEV_DIR}/RUST/rsoda` is a clean checkout at `ae76eee` (main) with push access; `gh` is authenticated.
- Python with scipy/radiopyo: `E:/miniconda3/envs/radiobio/python.exe`, run with `-I`.
