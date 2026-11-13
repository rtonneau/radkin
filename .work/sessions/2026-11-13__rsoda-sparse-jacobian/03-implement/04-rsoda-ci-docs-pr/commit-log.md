# Ticket 04: rsoda-ci-docs-pr

**Status:** ✅ Done

## Local Test Result

In `D:/DEV/RUST/rsoda`: `cargo fmt --check`, then for features `""` and `--features sparse`: `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release`, `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` → rc=0 for both (the same steps as the new CI). `gh pr view feat/sparse-jacobian --json url,headRefOid` → `https://github.com/rtonneau/rsoda/pull/3`, head `2f47645ff4aea19a53818befec4eeba6c0cc5f82`.

## Review Notes

Base: rsoda `81638e2`; commit `2f47645`.

- `.github/workflows/ci.yml` (new): Ubuntu and Windows × feature off and on; fmt, clippy `-D warnings`, tests, rustdoc `-D warnings`. Modelled on radkin's workflow.
- README: dependency note, a "Linear solvers and sparse Jacobians" section with an example (the same code runs as `tests/sparse.rs::readme_sparse_example`), and the test commands.
- Fixed a rustdoc broken link to `SparseLu` that appeared with the feature off.
- Version 0.2.0.
- Pushed the branch and opened PR rtonneau/rsoda#3 (linked in T3). Head SHA `2f47645` goes to ticket 05.

Review: inline.

## Blockers / Challenges

None. The CI results on GitHub weren't awaited here; the local run covered the same matrix steps on Windows.

## Commits

- ed436a2 chore(rsoda): ticket 04 CI, docs, 0.2.0 at rsoda 2f47645; PR rtonneau/rsoda#3

## Time Spent

2m (ticket-start.js to ticket-complete.js)
