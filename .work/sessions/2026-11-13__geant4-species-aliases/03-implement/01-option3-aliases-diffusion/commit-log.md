# Ticket 01: option3-aliases-diffusion

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`: all clean; lib tests 43 passed (including new `voxel::seed::tests::option3_names_all_covered`), all other suites ok.

## Review Notes

Checked names and D values against G4ChemDissociationChannels_option1.cc and the G4*.cc definitions; GetName forms match geant4-rs molecules.csv. 30 aliases (14 old + 16 new), 15 diffusion defaults (configuration-level D where option1 sets one, else definition D). Default ignore list stays empty, so water/bulk tags are not silently dropped. Existing test `parses_grid_and_diffusion` updated (OH_r default 2.8e-9 -> 2.2e-9). README diffusion/alias comments updated. Note: the coverage test lists 15 solute configurations + 5 non-solutes (H2O^0, H2O(B), H3Op(B), OHm(B), NoneM) = 20, not the ticket's 17; I followed the source.

Parent review (base 0365f79fc895db270b8bf2ac9470c956cce62f8f; pre-existing changes: gps session files only): independent reviewer skipped, localized low-risk change (two constant tables + one test) assessed directly. All 15 D values re-checked against G4ChemDissociationChannels_option1.cc and G4Hydrogen/G4Electron_aq/G4HO2/G4O2/G4O3/G4Oxygen.cc: match. Note: all 7 existing D values changed (not only OH_r). Parent fix: README alias example `"X^0" = "O2"` replaced by `"O2(B)" = "O2"`. Release tests re-run by parent: all pass. Accepted (implementer: sonnet, medium).

## Blockers / Challenges

None

## Commits

- 5574842 fix: alias every Geant4-DNA option3 species and use option3 diffusion defaults (ticket 01)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
