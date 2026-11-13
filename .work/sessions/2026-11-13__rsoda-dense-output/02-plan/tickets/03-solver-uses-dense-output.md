# Ticket 03: solver-uses-dense-output

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `Cargo.toml` rsoda `rev` points at the ticket-1 merge commit (git dependency, no `path`).
- [ ] `solve()` fills each `t_eval` point in `(t_prev, t_step]` from `nordsieck()` + the ticket-2 evaluator with `t_step` = the value returned by `step_tcrit`; `integrate(te)` no longer called; module doc updated.
- [ ] All radkin checks clean; 27-case acceptance matrix passes; `compare_radiopyo` on a `t_eval` case reported before/after.

**Files to Touch:**
- `Cargo.toml`, `Cargo.lock`
- `src/solver.rs`
- `README.md` / `docs` only if they describe the old t_eval path

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && pwsh tools/acceptance.ps1 -Python E:/miniconda3/envs/radiobio/python.exe -Radiopyo "$DEV_DIR/PYTHON/radiopyo"
```

Expected:
All checks pass; acceptance reports 27/27 cases within tolerance.

**Notes:**

Record the before/after max relative difference from `compare_radiopyo` on the labarbe 200-t_eval case in the commit log.
