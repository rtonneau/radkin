# Ticket 03: study-setup

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] An option3-matched radkin chemistry config: the reactions and rates of `G4EmDNAChemistry_option3.cc` read from the Geant4 sources, deaerated, with no background species. Every reaction is traced to its source line in comments.
- [ ] An aerated config: CONFIG_PAPER_01 (52 µM O2) plus the O(3P) reactions, as in `data/voxel_example.toml`.
- [ ] A converter script (Python `-I`) that writes radkin voxel TOMLs from the Geant4 outputs: event filter, a grid fitted to the seed bounding box plus a margin at a given voxel size, `energy_ev` from edep, and a chemistry choice.
- [ ] A benchmark of a short `radkin voxel` run per sweep point that predicts each full run's time and memory. The practical sweep points are recorded.

**Files to Touch:**
- `.scratch/tests/2026-11-13__geant4-single-track-gvalues/configs/*.toml`
- `.scratch/tests/2026-11-13__geant4-single-track-gvalues/tools/*.py`

**Verification Step:**

Run:
```bash
python -I .scratch/tests/2026-11-13__geant4-single-track-gvalues/tools/make_configs.py --check && cargo run --release --bin radkin -- voxel <a generated 2 um config> --out <scratch dir>
```

Expected:
The configs load and a coarse run finishes. Its `yields.csv` has G columns, with values at 1 ps that match the Geant4 seed G-values.

**Notes:**

The 1 ps G-values in radkin must equal the Geant4 seed counts divided by edep/100 eV. Use that to check the conversion.
