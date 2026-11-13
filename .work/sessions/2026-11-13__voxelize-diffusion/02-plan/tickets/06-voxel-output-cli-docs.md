# Ticket 06: voxel-output-cli-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Snapshots at requested times as .npy (one file per species per time) plus a JSON index.
- [ ] A CSV of domain-integrated amounts and G(t) (molecules/100 eV given the seed deposited energy if supplied).
- [ ] A `radkin voxel <config.toml>` CLI subcommand; an example config and a small example CSV in `data/`.
- [ ] `tools/compare_geant4.py` plots G(t) against a Geant4-DNA chem output; README section with voxel-size guidance and the mean-field caveat.
- [ ] fmt/clippy/test all clean.

**Files to Touch:**
- `src/voxel/output.rs`
- `src/bin/radkin.rs`
- `data/voxel_example.toml`
- `data/voxel_example_molecules.csv`
- `tools/compare_geant4.py`
- `README.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --bin radkin -- voxel data/voxel_example.toml
```

Expected:
All checks clean; the example run writes snapshots and a G(t) CSV.

**Notes:**

Energy deposited for G(t) comes from an optional seed field `energy_ev`.
