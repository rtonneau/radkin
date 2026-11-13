# Ticket 05: geant4-seed-import

**Status:** ✅ Done

## Local Test Result

`cargo test --release voxel` → 22 passed; `--test voxel_seed` → 1 passed; full suite ok; clippy clean.

## Review Notes

Base fbe988b. CSV reader (header-based, optional event_id, nm→m), alias table, ignore list, binning at 1e6/(N_A·V) µmol/L per molecule. Out-of-grid molecules are an error with the count and their bbox, or with `voxel.clip` they are dropped and the lost % is logged. `[[voxel.seed]]` (file relative to the config folder, time, events, offset_nm, energy_ev). `run_schedule` advances to each seed/output time, injects, and calls an output callback. Geant4 names verified in 11.4.3 source: `G4MolecularConfiguration` builds fName = definition name + '^' + charge (e.g. e_aq^-1, °OH^0, H3O^1, OH^-1, H_2^0); user IDs from `G4ChemDissociationChannels` (H3Op, OHm, ...). Both are mapped by default. Deviations: default aliases are a Rust const (`DEFAULT_GEANT4_ALIASES`) rather than `data/geant4_aliases.toml`, since overrides live in the run TOML anyway; non-dynamic species (e.g. H2O) must be ignored explicitly. Accepted inline.

## Blockers / Challenges

The schedule test first seeded e_aq in a 1 nm³ voxel (~1.7 M), where e_aq+e_aq→H2 polluted the H2 check: the mean-field caveat in action. Moved e_aq to an uninjected event.

## Commits

- 70af9ec feat(voxel): Geant4-DNA molecule CSV seeds, alias table, binning, clip policy and seed schedule (ticket 05)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
