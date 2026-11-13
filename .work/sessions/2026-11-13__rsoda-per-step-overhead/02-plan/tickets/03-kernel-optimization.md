# Ticket 03: kernel-optimization

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Based on the ticket-2 breakdown, `derive` and `jacobian`/`jacobian_sparse` are optimized where they matter (e.g. precomputed index layouts, fused rate and stoichiometry loops, avoiding redundant work per call), with the kernel staying allocation-free.
- [ ] A change is kept only if the median uninstrumented solve time drops by >= 5% on the cases it targets, with no case slower beyond noise. Otherwise it is reverted.
- [ ] All radkin golden tests (RHS/Jacobian/G-values at 1e-12, runs at rtol 1e-6) pass unchanged.
- [ ] docs/perf/README.md records each attempt, kept or dropped, with its before/after numbers.

**Files to Touch:**
- `src/kernel.rs`
- `src/matrices.rs`
- `docs/perf/README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && python -I tools/breakdown.py --repeat 10
```

Expected:
Checks are clean and all golden tests pass. The breakdown shows that the kernel share dropped for each kept change, or docs/perf records why nothing was kept.

**Notes:**

If the profile shows the kernel under ~15% of solve time on every case, record that and keep the ticket to a documented no-op. Floating-point reordering is allowed within the golden tolerances.
