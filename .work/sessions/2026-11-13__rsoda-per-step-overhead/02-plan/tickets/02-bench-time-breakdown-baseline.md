# Ticket 02: bench-time-breakdown-baseline

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] examples/bench.rs gains a `--breakdown` flag. It enables rsoda's `timing` feature through a radkin feature (e.g. `timing = ["rsoda/timing"]`, off by default) and wraps the kernel's rhs/jac/jac_sparse calls with timers. Per case, it prints kernel RHS, kernel Jacobian, factor, solve and remainder times, with counts.
- [ ] A script (`tools/breakdown.py`, or an extension of tools/benchmark.py that needs no radiopyo) runs the benchmark matrix, dense and sparse, and writes a markdown table.
- [ ] `docs/perf/README.md` holds the baseline breakdown per case, the uninstrumented median solve time, the per-RHS cost, the timer overhead, the machine/toolchain, and reproduction commands.
- [ ] The default build (no feature) is unchanged: fmt/clippy/tests clean, golden tests pass.

**Files to Touch:**
- `Cargo.toml`
- `examples/bench.rs`
- `src/solver.rs`
- `tools/breakdown.py`
- `docs/perf/README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo clippy --release --all-targets --features timing -- -D warnings && cargo test --release && python -I tools/breakdown.py --repeat 10
```

Expected:
Checks are clean, and the script prints a breakdown table covering every benchmark case, dense and sparse.

**Notes:**

While rsoda's ticket-1 commit is unpushed, use a local `[patch."https://github.com/rtonneau/rsoda"] rsoda = { path = "../rsoda" }` and do not commit it. Alternatively, push the rsoda branch and pin its rev. The rhs/jac timing wrapper can live in solver.rs behind `cfg(feature = "timing")`, or the bench can supply it. Keep the kernel allocation-free. Optimization decisions later use uninstrumented solve times.
