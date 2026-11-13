# Ticket 01: tonneau-config

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] data/tonneau_2025.toml loads, variable_GValues is true, and its beam is a single pulse (1e7 Gy/s, 30 Gy, LET 0.56).
- [ ] A test loads and runs it.
- [ ] README and data/README mention it.

**Files to Touch:**
- `data/tonneau_2025.toml`
- `tests/config.rs`
- `README.md`
- `data/README.md`

**Verification Step:**

Run:
```bash
cargo test --release --test config tonneau
```

Expected:
test result: ok

**Notes:**

Follow the paper_runs "original" scenario.
