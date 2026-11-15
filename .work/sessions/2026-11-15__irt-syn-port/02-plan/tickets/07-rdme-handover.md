# Ticket 07: rdme-handover

**Effort:** high

**Acceptance Criteria:**
- [ ] `radkin voxel` with an `[irt]` table and `solver = "ssa"` runs the IRT stage to the handover time, bins the molecules onto the voxel grid and continues with the RDME in the same process, one IRT trajectory per RDME trajectory.
- [ ] Molecule counts are conserved across the handover: the first RDME record equals the last IRT count for every species.
- [ ] The yields CSV covers both stages on one time axis.
- [ ] A seed scheduled after the handover time is a config error.
- [ ] Without `[irt]`, `radkin voxel` output is bit-identical to before.

**Files to Touch:**
- `src/irt/handover.rs`
- `src/irt/sim.rs`
- `src/ssa/voxel_sim.rs`
- `src/voxel/config.rs`
- `src/voxel/sim.rs`
- `src/bin/radkin.rs`
- `tests/irt_handover.rs`

**Verification Step:**

Run:
```bash
cargo test --release
```

Expected:
The whole test suite passes, including `irt_handover` and the existing `ssa_` tests unchanged.

**Notes:**

The deterministic Strang path is out of scope for the handover: `[irt]` with `solver = "lsoda"` in `radkin voxel` is refused with a clear message.
