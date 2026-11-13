# Ticket 01: rsoda-dky-nordsieck

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `Lsoda::dky(&self, t: f64, k: usize, out: &mut [f64]) -> Result<(), LsodaError>` evaluates DINTDY on the current history; returns an error before the first step, for `out.len() < n`, and for illegal `t`/`k` (`InputError::Interpolation`); never changes `t()`/`y()` or solver state.
- [ ] `Lsoda::nordsieck(&self) -> Nordsieck<'_>` (exported from the crate root) gives `yh`, `nyh`, `n`, `nq` (NQCUR), `nqu`, `tn`, `h` (HCUR), `hu`.
- [ ] Tests: `dky(tcur, 0)` equals `y()` after a step; `dky` at interior points matches `integrate(t)` bitwise; illegal t/k errors; `nordsieck()` fields agree with existing accessors.
- [ ] rsoda fmt/clippy/test clean; PR opened on rtonneau/rsoda and merged to main.

**Files to Touch:**
- `${DEV_DIR}/RUST/rsoda/src/lsoda.rs`
- `${DEV_DIR}/RUST/rsoda/src/lib.rs`
- `${DEV_DIR}/RUST/rsoda/src/interpolate.rs` (only if the view type lives there)
- `${DEV_DIR}/RUST/rsoda/README.md` (API mention)

**Verification Step:**

Run:
```bash
cd "$DEV_DIR/RUST/rsoda" && cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All checks pass; new dky/nordsieck tests listed as ok.

**Notes:**

Work in the rsoda repo on branch `feat/dense-output`, not in radkin. Merge commit SHA is needed by ticket 3.
