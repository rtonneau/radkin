# Ticket 05: walls-and-scheduled-tracks

**Effort:** xhigh

**Acceptance Criteria:**
- [ ] With a box, diffusion displacements and reaction and product positions are reflected at the walls.
- [ ] A pair within the cutoff radius of a wall also samples a reaction time against the partner's mirror image, across every wall within the cutoff radius, with the edge and corner images.
- [ ] Seeds carry a time and an offset, as the voxel solvers accept; a track is injected at its time and its molecules join the pair search then.
- [ ] A uniform distribution stays uniform between the walls (a chi-square test on slabs).
- [ ] The survival probability of one pair in a small box agrees with a brute-force Brownian walk with reflection within statistical error; the same test without the image rule fails, which shows the rule matters.
- [ ] Widely spaced tracks injected at scheduled times reproduce the single-track G(t) within statistical error.

**Files to Touch:**
- `src/irt/particles.rs`
- `src/irt/search.rs`
- `src/irt/solver.rs`
- `src/irt/sim.rs`
- `src/irt/config.rs`
- `tests/irt_walls.rs`
- `tests/irt_schedule.rs`

**Verification Step:**

Run:
```bash
cargo test --release irt_
```

Expected:
All `irt_` tests pass, including the three checks.

**Notes:**

The box is the voxel grid's box when the config has a `[voxel]` table, so that ticket 7 needs no second definition. G(t) of a run with scheduled tracks divides by the energy injected so far, as the RDME yields do.
