# Ticket 04: background-reactions

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_`: 29 passed, 0 failed (`irt_background` 5, `irt_model` 7, `irt_sample` 9, `irt_single_track` 8).

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure).

`python -I tools/compare_irt_geant4.py --events 100 --boundary outflow` with the switch off (the tool now writes `background_reactions = false`): 2281.1 reactions per run, RMS of z = 1.00, one value of the decade table outside 3 standard errors (H3O+ at 10 us). These are the numbers of ticket 03 to the digit (`.scratch/tests/2026-11-15__irt-syn-port/04-outflow/report.md`).

## Review Notes

Base: `162b31ddec73696622e5693dc8c4c7638fed430d`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- A reaction of one particle (first order, with a fixed species such as water, or with a background species) fires at an exponential time. One time is drawn per sampling pass for the sum of the rates of all particles, then the particle and its reaction in proportion to their rates. A pair reaction `A + B` whose partner has a background gets a twin of rate `k c` (and `2 k_obs c` for `A + A`, tested).
- The molecules of each background species consumed are counted (`StageStats::background_consumed`) and written to `irt.json` (`background_consumed`, with the background concentrations in `background`). No background is debited and no background particle exists.
- `[irt] background_reactions = false` turns off every reaction of one particle, first-order ones included. Default true. The comparison tool and the Geant4 tests set it to false.
- First-order decay: 4000 molecules follow `exp(-k t)` at 5 times within 4 standard errors. Scavenging of e_aq by 250 uM of O2, as an initial concentration and as a fixed concentration: same law.
- The Geant4 comparison with the switch off is unchanged to the digit: no random number is drawn when no such reaction exists.
- **Moved to ticket 05:** the uniform, well-mixed box against the unit-cell SSA or LSODA. It needs walls to keep the density uniform. A box large against the diffusion length is not an option with full resampling: every particle moves at every reaction, so the cost grows with the square of the particle count, and a box that stays uniform for a reaction half-life needs thousands of particles.

Other changes:
- The pair search now only sees particles of a species that has a reaction partner among the particles (inert products such as H2 are left out). Results are unchanged; the search of a run with one reactive species is no longer quadratic.
- The decay tests showed the cost of full resampling: every particle is moved at every event, so 20 000 independent molecules took 15 s. The tests use 2000 to 4000.

## Blockers / Challenges

None.

## Commits

- ab5985f feat: IRT single-molecule and background reactions (ticket 04)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
