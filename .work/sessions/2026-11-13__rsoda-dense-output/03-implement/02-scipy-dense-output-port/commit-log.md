# Ticket 02: scipy-dense-output-port

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → exit 0, 13 test binaries ok (log `.scratch/tests/2026-11-13__rsoda-dense-output/02-checks.log`). New `tests/dense.rs::dense_output_matches_scipy`: 361 real scipy LSODA steps (Robertson), 5 points each, 18 steps with a lowered order; worst relative difference 4.6e-15 (tolerance 1e-12 rel + 1e-15 of the step's max |y|).
Mutation check: disabling the `(HCUR/HU)^NQU` rescale fails the test at step 21 (3.5448575e-5 vs scipy 3.5448534e-5, ~1e-6 rel), so the golden exercises it and the scipy/DINTDY gap is truncation-level.

## Review Notes

Base f05b65eae3f1426e551fc2ac2f8d921814a0f501, clean tree apart from the gps session record. Acceptance: pure function `radkin::dense::lsoda_dense_output(&Nordsieck, t_step, t, out)`, no allocation, columns 0..=nqu, rescale of the last column when nq < nqu, `powf` for numpy `**` ✔. Golden: deviation from the plan, real scipy steps (raw `rwork[20:]`, `iwork[13:15]`, `rwork[10:12]` and `solver.dense_output()` values) instead of synthetic `yh` arrays. That covers `_dense_output_impl` end to end, including real order decreases and the final step where `t_step` (tcrit) != `tn`. `export_golden.py --only dense` regenerates just this file; the full export also writes it. Rust test at 1e-12 ✔.

## Blockers / Challenges

None. `tests/golden/dense_output.json` is 429 kB, the same order as runs.json.

## Commits

- a9aebb6 feat(dense): exact port of scipy LsodaDenseOutput, golden-tested (ticket 02)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
