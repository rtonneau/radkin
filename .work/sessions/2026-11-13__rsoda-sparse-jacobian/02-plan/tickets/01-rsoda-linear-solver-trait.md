# Ticket 01: rsoda-linear-solver-trait

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] rsoda exposes a public `LinearSolver` trait covering what DPRJA/DSOLSY do with the iteration matrix: form `P = I - h*el0*J` from a Jacobian source (user callback or FD of `rhs`), return `pdnorm` and the singular flag (`PrjaOutput`), and solve `P x = b` in place.
- [ ] Dense and banded (and the existing diagonal DSOLSY path) are implementations that keep the current arithmetic exactly (same loops, same `fnorm`/`bnorm`, same `dgefa`/`dgbfa`/`dgesl`/`dgbsl`).
- [ ] `stoda` and `lsoda` (`match_matrix_storage`, `set_jacobian`) use the trait; the `JacobianKind` to solver mapping lives in one place.
- [ ] Users can supply their own `LinearSolver` (e.g. an `OptionsBuilder`/`Lsoda` entry point that takes a boxed or generic solver), with a test using a custom solver that wraps the dense one.
- [ ] All existing rsoda unit, analytic and golden tests pass unchanged; no golden tolerance changes.

**Files to Touch:**
- `D:/DEV/RUST/rsoda/src/solsy.rs`
- `D:/DEV/RUST/rsoda/src/prja.rs`
- `D:/DEV/RUST/rsoda/src/stoda.rs`
- `D:/DEV/RUST/rsoda/src/lsoda.rs`
- `D:/DEV/RUST/rsoda/src/state.rs`
- `D:/DEV/RUST/rsoda/src/options.rs`
- `D:/DEV/RUST/rsoda/src/lib.rs`

**Verification Step:**

Run:
```bash
cd D:/DEV/RUST/rsoda && cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
fmt clean, no clippy warnings, every test passes (including `tests/golden.rs` and `tests/analytic.rs`) with no tolerance edits.

**Notes:**

Create branch `feat/sparse-jacobian` in the rsoda repo from `main` (312b742) before starting. Keep `PrjaInput`/`PrjaOutput` semantics (nje/nfe accounting, `hl0`, `y` restored after FD). The FD Jacobian logic needs `OdeSystem::rhs`, so the trait method takes the system as `&mut dyn OdeSystem` (or a generic) alongside `PrjaInput`. Also run radkin's tests against this rsoda through a temporary local `[patch]` (not committed) as an extra bit-identity check.
