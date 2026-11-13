# Ticket 06: pin-rsoda-and-final-report

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] The rsoda branch `perf/per-step-overhead` is pushed with its CI green and a PR open. rsoda's version is bumped (patch or minor) with a changelog/README note.
- [ ] radkin's Cargo.toml pins the pushed rsoda rev. No `[patch]` or `path` dependency remains.
- [ ] docs/perf/README.md holds the final breakdown and median solve times vs the baseline per case, plus the list of tried and dropped changes. The idea file `.work/idea/rsoda-per-step-overhead.md` is referenced.
- [ ] The radkin checks (fmt, clippy, tests) are clean.

**Files to Touch:**
- `Cargo.toml`
- `Cargo.lock`
- `docs/perf/README.md`
- `${DEV_DIR}/RUST/rsoda/Cargo.toml`
- `${DEV_DIR}/RUST/rsoda/README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && git grep -n "patch\|path = " Cargo.toml
```

Expected:
Checks are clean, and the grep shows no rsoda `[patch]` or `path` dependency.

**Notes:**

Pushing the rsoda branch and opening its PR is part of the approved design (rsoda changes are pushed to rsoda, then radkin bumps the rev). Do not merge the rsoda PR without the user.
