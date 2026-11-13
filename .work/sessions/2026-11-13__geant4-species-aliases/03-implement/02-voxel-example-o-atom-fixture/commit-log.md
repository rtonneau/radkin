# Ticket 02: voxel-example-o-atom-fixture

**Status:** ✅ Done

## Local Test Result

Command: `cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release && cargo run --release --bin radkin -- voxel data/voxel_example.toml --out .scratch/tests/2026-11-13__geant4-species-aliases/voxel_example` (exit 0).

- fmt and clippy -D warnings clean; 16 test binaries `ok`, 0 failed (golden tests included, unchanged); new `tests/voxel_geant4_seed.rs`: 3 passed.
- Voxel example: `radkin: voxel run of 768 voxels, 46 steps (0 rejected) in 0.6 s; 4 snapshots written`.
- Full log: `.scratch/tests/2026-11-13__geant4-species-aliases/evidence/02-verify.log`.

## Review Notes

- AC1: `data/voxel_example.toml` appends the 10 O(3P) k_reactions (values cross-checked against G4EmDNAChemistry_option3.cc lines 160-435; `O + O` written `2 O_atom -> O2`). Test `voxel_example_has_the_option3_o_atom_reactions` asserts 10 reactions with O_atom and that the H2O one gives 1.9e3 /s.
- H2O convention (radkin/radiopyo): H2O is a constant species ([H2O] = 55.933 mol/L at 25 degC from its polynomial) multiplied into the mass-action rate, so k is second order. `O_atom + H2O -> 2 OH_r` uses k_value = 33.97 (= 1.9e3 / 55.933), effective rate 1.9e3 /s within 1e-4; conversion stated in a TOML comment.
- Geant4 discrepancy noted in that comment: in G4EmDNAChemistry_option3.cc the comment says "1.90e3 / s" but the code sets `1.00e3 / s` for O + H2OB. Ticket value 1.9e3 used; parent may want to decide.
- AC2: CONFIG_PAPER_01.toml and golden data untouched; golden tests pass.
- AC3: fixture `tests/data/geant4_option3_event.csv` = event 1 of molecules.csv (9672 rows, 12 species). Event 1 is the only event holding all 12 species (HO_2^-1 occurs in 3 events, O^-1 in 10), so it is the smallest. Header `#` comment added, time_ps column kept. Test seeds it into the voxel_example chemistry with a 40^3 x 10 nm grid centred on the first `°O^0` via offset_nm (the track spans ~70 um), clip = true, `H2O^0` ignored: asserts placed + clipped + ignored = 9672, placed >= 300, clipped > 0, O_atom placed, and that the same seed with clip = false errors ("outside the grid").
- AC4/AC5 per the Verification Step above.

Parent review (base 55c990ff0b1bfc1609d3d7f9da6b7dcc7cdeb63e): independent reviewer skipped, data + test-only change assessed directly. Reactions and rates checked against the approved design; `2 O_atom -> O2` form is used by CONFIG_PAPER_01 too (`2 e_aq`). H2O convention (55.93 M multiplied in) documented and tested. Fixture 517 kB, provenance header present. New test re-run by parent: 3 passed. Accepted (implementer: sonnet, high).

## Blockers / Challenges

Positions in the fixture are tens of microns from the origin, hence the offset plus deliberate clipping in the test. The first draft of the H2O lookup in the reaction test matched H2O2 reactions; fixed by excluding H2O2. No review rounds yet.

## Commits

- b6e22cb feat: O_atom chemistry in voxel_example and a Geant4 option3 seed fixture test (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
