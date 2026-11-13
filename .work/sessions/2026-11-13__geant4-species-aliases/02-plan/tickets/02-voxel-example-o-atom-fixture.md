# Ticket 02: voxel-example-o-atom-fixture

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `data/voxel_example.toml` declares O_atom and the 10 Geant4 option3 O(3P) reactions: H_r+O_atom->OH_r 2.02e10; OH_r+O_atom->HO2 2.02e10; HO2+O_atom->O2+OH_r 2.02e10; O_atom+O_atom->O2 2.20e10; H2O2+O_atom->HO2+OH_r 1.60e9; H2+O_atom->H_r+OH_r 4.77e3; OH_minus+O_atom->HO2_minus 4.20e8; O2+O_atom->O3 4.00e9; HO2_minus+O_atom->O2_minus+OH_r 5.30e9 (M-1 s-1); O_atom+H2O->2 OH_r at 1.9e3 /s pseudo-first-order, written in radkin's H2O convention with a comment stating the conversion.
- [ ] radiopyo-derived configs (CONFIG_PAPER_01.toml etc.) are unchanged; golden tests unchanged and passing.
- [ ] A fixture CSV with one geant4-rs event that contains every species present in the run is committed under tests/ (smallest such event), and a test seeds it into the voxel_example chemistry (water ignored via policy) with zero errors and every row placed or deliberately clipped.
- [ ] `radkin voxel data/voxel_example.toml` still runs.
- [ ] fmt, clippy -D warnings and the release test suite pass.

**Files to Touch:**
- `data/voxel_example.toml`
- `tests/data/geant4_option3_event.csv` (new)
- `tests/voxel_geant4_seed.rs` (new, or an existing voxel test file)

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --bin radkin -- voxel data/voxel_example.toml --out .scratch/tests/2026-11-13__geant4-species-aliases/voxel_example
```

Expected:
All checks pass; the fixture test passes; the voxel example run finishes without error.

**Notes:**

Check radiopyo (${DEV_DIR}/PYTHON/radiopyo) for how H2O in a k_reaction is treated (constant [H2O] multiplied in or not) before writing O_atom+H2O. Fixture source: D:\DEV\GEANT4\SIM\geant4-rs\out\molecules.csv; positions are in nm, possibly far from the origin: use a grid/offset in the test that covers the event (or clip), keeping the test fast. A config including CONFIG_PAPER_01 may need extra reactions appended after the include; check the merge rules in src/config/.
