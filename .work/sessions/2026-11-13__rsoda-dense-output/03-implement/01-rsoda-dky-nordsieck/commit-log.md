# Ticket 01: rsoda-dky-nordsieck

**Status:** ✅ Done

## Local Test Result

rsoda (`feat/dense-output`): `cargo test --release` → lib 103 passed (5 new: dky_before_the_first_step_is_an_error, dky_at_the_mesh_point_is_the_step_output, dky_inside_the_last_step_matches_integrate_and_changes_nothing, dky_rejects_illegal_t_k_and_short_output, nordsieck_view_agrees_with_the_accessors), integration 10+1+1 passed, doctests ok; `cargo doc --no-deps` 0 warnings. `cargo clippy -D warnings`: 17 errors, identical count on `main` (all in linalg/norms/options/prja, untouched). `cargo fmt --check` fails on `main` too (rsoda has never been rustfmt-formatted), so no fmt run on the PR.
radkin after the rev bump: `cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release` → exit 0 (log `.scratch/tests/2026-11-13__rsoda-dense-output/01-radkin-checks.log`).

## Review Notes

Base radkin HEAD 4384c71c941e33fc29e3adfac29d3887dabddc1b; rsoda base ae76eee69d6d494be9253eaee819c77ee962c934, clean. Acceptance criteria: dky signature/errors/no state change ✔ (tests), `Nordsieck` exported from crate root with yh/nyh/n/nq/nqu/tn/h/hu ✔, tests ✔, PR https://github.com/rtonneau/rsoda/pull/2 merged → merge commit 312b742.
Scope change: the radkin rev bump (Cargo.toml, Cargo.lock, `RSODA_REV` in src/output.rs) moved from ticket 3 into this ticket, because ticket-complete only records radkin commits; the API is additive, so radkin stays green. Ticket 3 keeps the solver change.

## Blockers / Challenges

A first `cargo fmt` in rsoda reformatted all 16 files (no rustfmt config, code formatted wider than the default); reverted with `git checkout -- .` and the change re-applied by script without fmt. Pre-existing rsoda clippy/fmt debt left as is (out of scope; touching linalg could change a faithful port).

## Commits

- 3f225bb chore(deps): bump rsoda to 312b742 (dky and Nordsieck view) (ticket 01)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
