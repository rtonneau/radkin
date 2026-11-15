# Ticket 03: single-track-solver

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_`: 24 passed, 0 failed (`irt_model` 7, `irt_sample` 9, `irt_single_track` 8). `cargo test --release --lib irt`: 8 passed (config, search, species).

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure).

`python -I tools/compare_irt_geant4.py --events 100` (free space), exit code 1:
- RMS of z = 1.76 over 1065 values; 6 values of the decade table outside 3 standard errors, all at 1 us or 10 us (e_aq, OH, H, H3O+, OH- at 10 us; H3O+ at 1 us).
- With `--t-max 3.2e-7`: RMS of z = 0.93 over 825 values, every value of the decade table within 3 standard errors, exit code 0.

`python -I tools/compare_irt_geant4.py --events 100 --boundary outflow` (Geant4's 2 um box), exit code 1:
- RMS of z = 1.00 over 1065 values; every value of the decade table from 10 ps to 10 us within 3 standard errors except H3O+ at the last record (-3.2).

Report: `.scratch/tests/2026-11-15__irt-syn-port/03-report.md`, with the tool's own reports in `03-free/`, `03-free-to-300ns/`, `03-outflow/` and the checks in `03-check-*/`.

## Review Notes

Base: `7600697cbd851b716e9b3b2cb2d65cf0e85f3e4c`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- `src/irt/` holds the particles and the event loop with full resampling (`solver.rs`), the pair search on a grid of cells with a cutoff radius that shrinks with the time left (`search.rs`), the `[irt]` config (`config.rs`), the run and its ensemble (`sim.rs`) and the outputs (`output.rs`).
- `radkin irt <config> --out <dir>` writes `yields.csv` (the columns of the RDME yields: `t`, `N_`, `G_`, `N_..._std`, `G_..._std`) and `irt.json`. One trajectory per random stream, in parallel with rayon.
- A config with a `[beam]` and `[irt]` is refused.
- Counts change only through reactions: a run of inert molecules keeps every count; on the Geant4 network the change of every count equals the sum of the reactions that fired, and the charge is constant at every record.
- A regression test runs the ten reference events of `tests/data/irt/` (32 trajectories each) and asserts every species within 4 standard errors of Geant4 up to 1 us; largest |z| 1.19.
- **Acceptance against Geant4: met with a qualification.** The criterion asked for free space from 1 ps to 10 us. The Geant4 reference is not free space: it ran in a 2 um water box and loses molecules at its faces from about 300 ns on (its net charge drifts from -5 to -83, its H2 decreases). In free space radkin agrees to 300 ns (RMS of z 0.93). With the same box (`boundary = "outflow"`) it agrees to 10 us (RMS of z 1.00), except H3O+ in the last half decade, where the result depends on Geant4's transport steps, which radkin does not port. A reference whose molecules cannot leave the water would settle the last decade in free space.

Choices made while implementing (all stated in the code docs):
- `boundary = "outflow"` was added to reproduce Geant4's box: a molecule that a move puts outside the box is removed. It is not an absorbing wall. The plan's mitigation named a "position clamp"; the data show Geant4 removes the molecules.
- `start_contact` (default true): a pair injected inside its effective radius reacts at once with the contact probability. It must be false for seeds Geant4 wrote after its chemistry started (1.001 ps); true on these seeds gives +3.1 standard errors of H2O2 at 10 ps.
- A pair found inside its effective radius after a move gets no reaction time in that pass (as planned), and bounds the next step to `r_eff^2 / D`, after which every particle is moved and looked at again. Without that bound such a pair would wait for a reaction elsewhere. A shorter bound (the depth of the pair inside its radius) made the run crawl; treating the pair as in contact disagrees with Geant4 at 10 ps (RMS 1.50).
- The random seed key is `random_seed`, since `[[irt.seed]]` already names the molecule seeds.
- The record times are not stops: particles are moved at reactions (and at the end), and a record reads the counts.
- `[irt]` reuses the seed, alias and box definitions of `[voxel]` when it has none of its own; `seed_spec` and the Geant4 alias parsing of `src/voxel/config.rs` are now shared helpers (`radkin voxel` behaviour unchanged, its tests pass).

Not in the ticket's file list: `src/voxel/config.rs` (shared helpers), `src/irt/model.rs` (reaction count, extra diffusion map), `tests/data/irt/README.md`. `src/irt/particles.rs` is not a file of its own: the particle state is part of `solver.rs`. `src/error.rs` did not need a change.

## Blockers / Challenges

- The Geant4 reference box (above). Three rules for the faces were measured before keeping the simplest one.
- The seed pass and the G(t) pass of Geant4 are two realisations of the pre-chemical stage; the first version of the comparison ignored it and reported H atoms 10 standard errors off at 1.26 ps. The tool now adds the difference of the initial states to the error.
- Geant4's time stepper does not clear the reaction set after each reaction step, only after a step limited by a process (`Prepare` tests `IsInteractionStep`). Its IRT_syn is therefore closer to incremental resampling in an unbounded medium, and resamples everything only when a molecule near a face limits the step. radkin's full resampling still matches it statistically. This matters for ticket 06.

## Commits

- 19d52db feat: IRT single-track solver with full resampling, radkin irt (ticket 03)

## Time Spent

29m (ticket-start.js to ticket-complete.js)
