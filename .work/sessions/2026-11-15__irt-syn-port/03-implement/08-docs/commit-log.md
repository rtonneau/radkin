# Ticket 08: docs

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`: the three checks are clean (33 test binaries `ok`, no failure).

The README example runs: `radkin irt data/irt_example.toml --out <dir>` prints "irt run of 66 molecules, 32 trajectories, 762 pair, 24 contact and 12 single-molecule reactions".

## Review Notes

Base: `a305ab0607e6c486531e233f9c7178971a8d6460`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- ADR 0007 (`.work/adr/0007-irt-stage.md`) records the IRT stage: own solver with Geant4's formulas, the two resampling rules, the type rule, the background, the image rule at walls and its limit, one handover time, and what is not done.
- `docs/irt.md` describes the model, the kernel, the two pair rules, background reactions, the three boundaries, seeds and scheduled tracks, the handover, the `[irt]` config, the outputs, the validation against Geant4, the cost and the limits (independent pairs, provisional multi-track and aerated results).
- `README.md` has an "IRT stage" section; `AGENTS.md` has the module in its layout, the `irt_` test prefix, the tools and the reference data; `tools/README.md` lists `export_g4_chem.py` and `compare_irt_geant4.py`.
- The glossary edit of the grill was committed with the plan (`aee437c`). This ticket corrects one of its entries (full resampling is not what Geant4's IRT_syn does) and adds the terms settled while building: pair reaction, molecule reaction, background species, contact rule, image pair, outflow box.

Numbers in the docs come from the saved reports: `.scratch/tests/2026-11-15__irt-syn-port/03-report.md` (Geant4 comparison), the tool reports of tickets 03 and 06, the test output of ticket 05 (walls) and `docs/perf/irt-scaling.md`.

## Blockers / Challenges

None.

## Commits

- 404834b docs: IRT stage (docs/irt.md, README, AGENTS, tools README) (ticket 08)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
