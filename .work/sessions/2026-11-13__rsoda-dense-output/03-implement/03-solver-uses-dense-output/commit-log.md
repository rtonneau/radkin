# Ticket 03: solver-uses-dense-output

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → exit 0, 13 test binaries ok (`03-checks.log`). `pwsh tools/acceptance.ps1 -Python E:/miniconda3/envs/radiobio/python.exe -Radiopyo $DEV_DIR/PYTHON/radiopyo` → 27/27 cases passed (`03-acceptance.log`).
compare_radiopyo before/after (`03-compare-before.txt`, `03-compare-after.txt`), Labarbe run with 200 t_eval points: no beam, worst species rel. diff 2.1e-08 → 2.1e-08; constant beam, 2.5e-08 → 2.5e-08 (rates 1.0e-08 → 1.0e-08). The radkin output did change: 198/200 points differ, up to 7.4e-08 relative (DINTDY vs scipy formula). At the 2-digit precision the tool reports, the gap to radiopyo is the same, so stepping differences between rsoda and scipy's Fortran LSODA dominate what's left.

## Review Notes

Base d7edb9c14cf4e99a45abda5e2030689b3a47c8bd. Acceptance: the rev bump was done in ticket 01 (312b742; see that log) ✔. `solve()` fills every t_eval point of a step from `nordsieck()` + `lsoda_dense_output` with `t_step` = the value `step_tcrit` returned, including a point equal to the step end (scipy evaluates `sol(t)` there too). `integrate(te)` is gone and the module doc is updated ✔. Checks, acceptance and compare ✔. `docs/superpowers/plans/2026-11-13-radiopyo-port.md:195` still describes the old path; it is the historical Phase-1 plan, so I left it.

## Blockers / Challenges

None.

## Commits

- 744b1de feat(solver): fill t_eval with scipy dense output, drop integrate(te) (ticket 03)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
