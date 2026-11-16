# Ticket 08: docs

**Effort:** medium

**Acceptance Criteria:**
- [ ] ADR 0007 records the IRT stage: own solver with Geant4's sampling formulas, the two resampling rules, the type rule, the image rule at walls, one handover time.
- [ ] `docs/irt.md` describes the algorithm, the `[irt]` config, the outputs, the acceptance result against Geant4, the measured scaling, and the provisional status of multi-track aerated results.
- [ ] `README.md`, `AGENTS.md` (layout, `irt_` test prefix) and `tools/README.md` mention the IRT stage and its tools.
- [ ] The glossary edit of the grill is committed.

**Files to Touch:**
- `.work/adr/0007-irt-stage.md`
- `.work/GLOSSARY.md`
- `docs/irt.md`
- `README.md`
- `AGENTS.md`
- `tools/README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
The three checks are clean.

**Notes:**

Numbers in the docs come from the reports saved in the session scratch directory, not from memory.
