# Tools (Python, radiopyo conda env)

- `export_golden.py`: converts radiopyo's golden data (`tests/golden/data/*` of the radiopyo
  checkout) plus a few reference values computed with radiopyo into `tests/golden/*.json`.
- `compare_radiopyo.py`: runs radiopyo and radkin on the same inputs and compares the outputs.
- `acceptance.ps1`: the Phase-1 acceptance matrix built on `compare_radiopyo.py`.
- `benchmark.py`: runtime comparison on 9 cases. Both codes are timed in-process
  (load = build from TOML, solve = integration, no output files); radkin side is
  `examples/bench.rs`. `python -I tools/benchmark.py --repeat 5 --out bench.json`.
- `export_g4_chem.py`: turns the reaction table Geant4-DNA dumped at run time
  (`tests/data/irt/g4_chemistry_dump.txt`) into `data/chem_g4.toml` and
  `data/chem_g4.diffusion.json`, with the Geant4 reaction type of every pair reaction
  as `irt_type`. `python -I tools/export_g4_chem.py [DUMP [OUT.toml]]`. No radiopyo
  needed.
- `compare_irt_geant4.py`: runs radkin's IRT stage (`radkin irt`) on every event of a
  Geant4-DNA IRT_syn reference and compares G(t) of every species at Geant4's checkpoint
  times, with the standard error of the difference. `python -I
  tools/compare_irt_geant4.py --g4 <run dir> --same-pass --events 2000` for a run whose
  seeds were written by the G(t) pass itself (`/rs/chem/seedDump true`), or `--events 100
  --boundary outflow` for the first reference (seeds of a separate pass, in a 2 um box).
  It needs a release build of radkin and the Geant4 reference directory; no radiopyo.
  See `docs/irt.md`.
