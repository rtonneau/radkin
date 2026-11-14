# Ticket 02: chemistry-export-and-tooling

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] A script converts the Geant4 dump into a radkin TOML (species, D, reactions, rates), deaerated, with each mapping commented
- [ ] An aerated config: CONFIG_PAPER_01 with 52 µM O₂ plus the O(³P) reactions
- [ ] Seed conversion and run scripts build per-event grids (bounding box + 2 µm) starting at the handover time and running to 10 µs
- [ ] A plot script pools events 0–4 and overlays the old IRT G(t)
- [ ] A 1 µm smoke run per chemistry runs, and G at the handover time matches the seed G

**Files to Touch:**
- `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/tools/*`
- `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/configs/*`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/tools/run_sweep.py --only smoke
```

Expected:
Exit 0 for each smoke case. G at the handover time equals the seed G.

**Notes:**

Reuse the previous session's tools where possible. Commit radkin code only if a bug is found, with a test.
