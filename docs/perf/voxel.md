# radkin performance: voxel runs (Phase 2)

Benchmark and accuracy gate for the voxel reaction-diffusion run (gps session
`2026-11-14__reduce-computational-burden`). Later changes that cut the cost of a
voxel run are measured with it: a lever stays only if it is faster and the
accuracy columns below do not get worse beyond what the change is meant to trade.

## How it is measured

`examples/voxel_bench.rs` runs the Geant4-seeded case `data/voxel_example.toml`
(8x8x12 voxels of 6 nm, seeds `data/voxel_example_molecules.csv`, CONFIG_PAPER_01
chemistry plus O(3P)) `--repeat` times in-process. Config loading is not timed.
The end time is the configuration's `t_end = 1e-7 s` (override with `--t-end`):
about 0.6 s per repeat on 16 threads, 3.4 s on one thread. A snapshot is always
taken at `t_end`; the final state is what is compared.

- **Counters** come from `radkin::voxel::reaction::stats()` (read-only atomics in
  the reaction step, no effect on the integration): voxel half-steps (all voxels,
  rejected attempts included), background voxels among them, rsoda solves actually
  run (active voxels plus the shared background integration), and their summed
  `nfev` and `nst`. Accepted/rejected splitting steps come from `VoxelOutput`.
  `njev` is not reported (it read 0 in these solves, not verified further).
- **Reference**: the same case with the splitting control tightened to
  `rtol = 1e-6`, `atol = 1e-6` umol/L (the case uses `rtol = 1e-2`, `atol = 1e-3`;
  atol is tightened too because at 1e-3 it would dominate the small species).
  Per-voxel rsoda tolerances are unchanged (`1e-10` / `1e-12`). It is written by
  `--write-reference` to
  `.scratch/tests/2026-11-14__reduce-computational-burden/voxel_reference.json`
  (not committed; 578 steps, about 5 s) and checked against the case (species,
  grid, `t_end`) before use.
- **Accuracy**: for each species the relative error of the domain total (sum over
  voxels of the final state; denominator floored at the run's `atol`), and the max
  over all voxels and species of `|y - y_ref| / (atol + rtol * |y_ref|)` with the
  run's own control tolerances (below 1 = within what the run asked for).
- The last stdout line is `VOXEL_BENCH key=value ...` for scripts.

Reproduce:

```bash
cargo run --release --example voxel_bench -- --write-reference   # once per case / t_end
cargo run --release --example voxel_bench -- --repeat 3
RAYON_NUM_THREADS=1 cargo run --release --example voxel_bench -- --repeat 3
```

## Baseline

radkin `1465054`, rsoda `f5d1312`. Intel Core i7-7820X 3.6 GHz (8 cores, 16
threads), Windows 11, rustc 1.99.1, `--repeat 3`, `parallel` feature (default).

| threads | median [s] | min [s] |
|---|---|---|
| 16 (default) | 0.574 | 0.570 |
| 1 (`RAYON_NUM_THREADS=1`) | 3.365 | 3.364 |

Work (identical at any thread count, deterministic):

| quantity | value |
|---|---|
| splitting steps accepted / rejected | 46 / 0 |
| voxel half-steps (of which background) | 211968 (29039) |
| rsoda solves | 183205 |
| nfev (sum over solves) | 4389691 |
| nst (sum over solves) | 2202873 |

About 24 rsoda RHS calls per solve and 0.77 us of wall time per RHS call per
thread (3.365 s / 4.39M).

Accuracy of the final state against the reference (case rtol 1e-2, atol 1e-3):

| metric | value |
|---|---|
| max per-voxel scaled error | 0.193 |
| max domain-total relative error | 1.50e-3 (HO3) |
| other domain totals | 4.9e-5 (H2O2) to 1.1e-3 (O3); e_aq 5.5e-4, OH_r 3.4e-4 |

The reference itself is the 5.2 s run (578 steps), so the default case is about 9x
cheaper than the reference while staying inside its requested tolerance.
