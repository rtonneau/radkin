# Ticket 05: geant4-seed-import

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] A CSV reader for molecules (species, x, y, z in nm, optional event_id); header-based columns.
- [ ] A `[geant4.aliases]` table plus a shipped default for the standard Geant4-DNA water species; unknown names are an error unless listed in `ignore`.
- [ ] Binning at 1/(N_A·V) per molecule; out-of-grid molecules are an error with count and bbox, or are clipped with the lost fraction reported.
- [ ] A `[[seed]]` schedule (file, event filter, time, offset), injected by the driver at exact times.
- [ ] Tests: the binned total molecule count is preserved; alias errors; clip reporting; two seeds injected at different times.

**Files to Touch:**
- `src/voxel/seed.rs`
- `src/config/mod.rs`
- `src/voxel/driver.rs`
- `data/geant4_aliases.toml`
- `tests/voxel_seed.rs`

**Verification Step:**

Run:
```bash
cargo test --release seed
```

Expected:
Seed tests pass.

**Notes:**

Verify Geant4-DNA molecule names against the G4 source (G4MoleculeTable definitions), not from memory.
