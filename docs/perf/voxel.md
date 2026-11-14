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
  `ReactionStats::njev` (added with ticket 02, not printed by the bench) is 0
  over all solves: rsoda counts a Jacobian evaluation in every BDF matrix
  setup, so no voxel solve of this case ever leaves the Adams (non-stiff)
  method, and the analytic Jacobian and dense LU are never used.
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

## Lever 1: lean voxel solve (ticket 02) - kept

The reaction half-step integrates each voxel with `solver::integrate_to`
instead of `solver::solve` with `t_eval = [t + h]`:

- no output storage (`Solution`, `t_eval` and options clones), the result is
  written in place;
- the kernel `Workspace` (and sparse pattern) is built once per rayon pool
  thread per half-step and reused for all voxels of that thread (one slot per
  `rayon::current_thread_index()`). rayon's `try_for_each_init` was tried
  first: it builds one scratch per split job, which was 154224 builds for
  183205 solves at 16 threads (2316 at 1 thread), so it gave nothing in
  parallel;
- each voxel starts with `CARRY_FACTOR = 0.01` times the last step its previous
  half-step took (not the step shortened to land on the end), carried along the
  accepted trajectory in `VoxelRun` (reset on injection, not in checkpoints).

`integrate_to` with the same first step is bit-identical to `solve` (unit test),
and results do not depend on the thread count (tests with 1- and 4-thread
pools). rsoda restarts at order 1 on every solve, so the carried step mostly
saves failed or tiny first steps: the factor scan gave nfev 4.65M (1.0),
4.34M (0.1), 4.27M (0.03), 4.22M (0.01), 4.39M (0.003); rsoda's own initial
step (`h0 = 0`) gave 4.63M. 0.01 is used.

A/B on the same machine, alternating the ticket-01 build (`99e8d37`) and this
one; medians of the per-run medians (16 threads: 6 runs of `--repeat 5`;
1 thread: 3 runs of `--repeat 3`). The ticket-01 build ran a bit faster at 16
threads than in the baseline table above (0.541 vs 0.574 s).

| threads | before [s] | after [s] | change |
|---|---|---|---|
| 16 | 0.541 | 0.505 | -6.7% |
| 1 | 3.410 | 2.885 | -15.4% |

| quantity | before | after |
|---|---|---|
| rsoda solves | 183205 | 183205 |
| nfev | 4389691 | 4219543 |
| nst | 2202873 | 2174060 |
| njev | 0 | 0 |
| splitting steps accepted / rejected | 46 / 0 | 46 / 0 |
| max per-voxel scaled error | 0.193 | 0.193 |
| max domain-total relative error | 1.496e-3 (HO3) | 1.496e-3 (HO3) |

The other domain totals are unchanged at the printed precision (HO2 6.106e-5
-> 6.107e-5).

Not done, needs rsoda: reusing rsoda's own buffers (the `State` arrays and the
dense `n x n` iteration matrix, allocated by every `Lsoda::new` and never used
here since njev is 0). A local prototype (`Lsoda::restart(t0, y0, h0)` that
resets the state in place and keeps the buffers and linear solver, used from
the same per-thread scratch; not committed) was bit-identical in the unit tests
and measured, in a separate A/B session where the ticket-01 build ran at
0.564 s (16 threads) / 3.71 s (1 thread) and this lever at about 0.53 s /
3.14 s, 0.495 s at 16 threads (about -7% on top of this lever) and 2.94 s at
1 thread (no gain there). Left as a follow-up idea
(`.work/idea/rsoda-lsoda-restart.md`).

## Lever 2: fused junction half-reactions (ticket 03) - kept

In the dt/2 path of `VoxelRun::advance_to`, the two half-step Strang steps are
run by `strang_pair`: reaction dt/4, diffusion dt/2 (Rannacher-damped only when
the attempt is), reaction dt/2 (the second half-reaction of the first sub-step
and the first of the second, fused), diffusion dt/2 (never damped), reaction
dt/4. That is 5 voxel reaction solves per attempted step instead of 6. The
background integration and the per-voxel carried steps follow the same calls
(the carry is committed only on accepted steps, as before). Diffusion is
unchanged; the fused reaction replaces two consecutive solves with no operator
in between, so the splitting error does not grow (the gate is unchanged).

A/B, alternating the ticket-02 build and this one (16 threads: 4 runs of
`--repeat 5`; 1 thread: 3 runs of `--repeat 3`), medians of per-run medians:

| threads | before [s] | after [s] | change |
|---|---|---|---|
| 16 | 0.512 | 0.434 | -15.1% |
| 1 | 2.890 | 2.518 | -12.9% |

| quantity | before | after |
|---|---|---|
| voxel half-steps (background) | 211968 (29039) | 176640 (24360) |
| rsoda solves | 183205 | 152510 |
| nfev | 4219543 | 3693056 |
| nst | 2174060 | 1896047 |
| splitting steps accepted / rejected | 46 / 0 | 46 / 0 |
| max per-voxel scaled error | 0.193 | 0.193 |
| max domain-total relative error | 1.4963e-3 (HO3) | 1.4963e-3 (HO3) |
