# Ticket 08: rdme-tau-leaping

**Acceptance Criteria:**
- [ ] `tau_leap = true` with `radkin voxel` leaps reactions and hops per voxel with Poisson counts; a channel is critical (fired exactly, one at a time) when the population it consumes is below 10
- [ ] The step is chosen with the Cao–Gillespie–Petzold bound over every (voxel, species) population, hops in and out included; when it is below 10/a0 the run falls back to exact next-subvolume events
- [ ] Counts never go negative; a rejected leap is halved and redone; molecule numbers are conserved by hops
- [ ] A leap never crosses a seed time or a snapshot time; coarsening still works
- [ ] Fixed-seed test: tau-leaped and exact RDME ensembles agree on species totals within 4 combined SEM on a large-count grid; speed-up recorded in `docs/ssa.md`
- [ ] fmt, clippy and tests pass

**Files to Touch:**
- `src/ssa/rdme.rs`
- `src/ssa/rdme_tau.rs`
- `src/ssa/voxel_sim.rs`
- `src/ssa/mod.rs`
- `tests/ssa_rdme.rs`
- `docs/ssa.md`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` tests pass, including the voxel tau-leaping test.

**Notes:**

A hop is a first-order channel `X(v) -> X(w)` of propensity `n·D/h²` per face. After a leap the next-subvolume queue is rebuilt before the next exact phase.
