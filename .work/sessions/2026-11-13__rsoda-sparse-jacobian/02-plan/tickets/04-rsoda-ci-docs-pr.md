# Ticket 04: rsoda-ci-docs-pr

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] `.github/workflows/ci.yml` runs fmt, clippy (`-D warnings`) and tests on Ubuntu and Windows, with the `sparse` feature off and on.
- [ ] README documents the `LinearSolver` trait, the sparse Jacobian kinds and the `sparse` feature, with a short example; rustdoc builds without warnings.
- [ ] The crate version is bumped to 0.2.0 (a breaking API change).
- [ ] Branch `feat/sparse-jacobian` is pushed to rtonneau/rsoda and a PR to `main` is opened; its head commit SHA is recorded for ticket 05.

**Files to Touch:**
- `D:/DEV/RUST/rsoda/.github/workflows/ci.yml` (new)
- `D:/DEV/RUST/rsoda/README.md`
- `D:/DEV/RUST/rsoda/Cargo.toml`

**Verification Step:**

Run:
```bash
cd D:/DEV/RUST/rsoda && cargo doc --no-deps --features sparse && git status --short && gh pr view feat/sparse-jacobian --json url,headRefOid
```

Expected:
docs build without warnings, a clean tree, and the PR URL with its head SHA.

**Notes:**

Mirror radkin's `.github/workflows/ci.yml` layout. Link the rsoda PR in the host after creating it.
