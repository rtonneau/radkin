# Ticket 02: configs-and-pilot

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Config generator writes, per (chemistry, handover, voxel): `solver = "ssa"`, default mesh correction, 4 trajectories, grid equal to the 2×2×2 µm water box, the seed at its handover time, and outputs to 10 µs. The chemistries are chem_g4 deaerated, chem_g4 + 52 µM O2, and CONFIG_PAPER_01 aerated with O(³P).
- [ ] Pilot ran: event 0, 1 ps, 10 nm, deaerated. Wall time, peak memory and event counts are recorded in `.scratch/.../pilot.md`.
- [ ] Aerated O2 species level chosen from radkin's suggestion and checked against level 0 on a short pilot. The choice is recorded.
- [ ] The sweep plan is sized from the pilot: trajectories in parallel, plus any reduction for 5 nm.

**Files to Touch:**
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/make_configs.py`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/configs/`
- `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/pilot.md`

**Verification Step:**

Run:
```bash
cargo build --release && python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/make_configs.py --check
```

Expected:
All configs load (radkin dry check or a short run). The pilot's yields.csv covers 1 ps to 10 µs, and its G at 1 ps equals the seed G.

**Notes:**

Reuse chem_g4.toml from `.scratch/tests/2026-11-14__geant4-radkin-irt-syn-gvalues/configs/`, and the seed CSV format and aliases from earlier sessions.
