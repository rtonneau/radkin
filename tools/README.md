# Tools (Python, radiopyo conda env)

- `export_golden.py`: converts radiopyo's golden data (`tests/golden/data/*` of the radiopyo
  checkout) plus a few reference values computed with radiopyo into `tests/golden/*.json`.
- `compare_radiopyo.py`: runs radiopyo and radkin on the same inputs and compares the outputs.
- `acceptance.ps1`: the Phase-1 acceptance matrix built on `compare_radiopyo.py`.
