# Ticket 03: report-and-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] The single-track 10 us G-value comparison is rerun with the fix. The report under `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/` compares RDME with the deterministic voxel run and the earlier numbers.
- [ ] `docs/ssa.md` describes the reversible-pair correction and the background run, including its cost and the opt-out flag.
- [ ] ADR 0006 is consistent with what was implemented.

**Files to Touch:**
- `docs/ssa.md`
- `.work/adr/0006-rdme-reversible-pair-correction.md`
- `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/report.md`

**Verification Step:**

Run:
```bash
cargo test --release && git diff --stat
```

Expected:
Tests pass. The docs and the report have been written.

**Notes:**

Reuse the inputs and scripts of `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/`. If the full 5 nm run is too heavy, run 20 and 10 nm and say so in the report.
