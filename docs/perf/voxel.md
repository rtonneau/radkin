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

## Lever 3: embedded splitting error estimate (ticket 04) - dropped

Goal: replace step-doubling (5 reaction solves and 3 diffusions per attempt
since lever 2) with an error estimate built from fewer solves of the same
step. Candidates, per attempted step of `dt`:

- **Strang vs Lie**, sharing the first half-reaction: `y1 = R(dt/2) y0`, then
  Strang `D(dt), R(dt/2)` and Lie `R(dt/2), D(dt)`. 3 solves of dt/2, 2
  diffusions. The estimate is the Lie error, O(dt^2): controller exponent -1/2.
  Strang is kept.
- **Strang vs reversed Strang**: `R/2 D R/2` (S1) and `D/2 R D/2` (S2), error
  `K * |S1 - S2|`, exponent -1/3. 3 solves (two of dt/2, one of dt) and 3
  diffusions. S1, S2 or their average is kept. One solve over dt costs about
  1.22 solves over dt/2 (S1 2.07M nfev, S2 1.26M over the 60 steps), so this
  is about 1.6 Strang steps of work per attempt.

The candidates ran from one exploration build (estimator chosen by an
environment variable; not committed). The work counts are deterministic.
Official case (t_end 1e-7, gate 1.4963e-3 / 0.193):

| estimator, kept solution | K | steps | rsoda solves | nfev | max total err | max scaled err | gate |
|---|---|---|---|---|---|---|---|
| step-doubling (lever 2) | - | 46 | 152510 | 3693056 | 1.496e-3 | 0.193 | - |
| Strang vs Lie, Strang | - | 135 | 293772 | 5425416 | 5.98e-4 | 0.133 | pass |
| reversed, S1 | 1 | 60 | 123903 | 3358046 | 2.656e-3 | 0.348 | fail |
| reversed, average | 1 | 60 | 123910 | 3362639 | 1.147e-3 | 0.176 | pass |
| reversed, average | 2/3 | 54 | 110077 | 3069152 | 1.481e-3 | 0.229 | fail |
| reversed, S2 | 1 | 60 | 123909 | 3339006 | 9.37e-4 | 0.120 | pass |
| reversed, S2 | 2/3 | 54 | 110074 | 3049868 | 1.210e-3 | 0.154 | pass |
| reversed, S2 | 1/2 | 51 | 103156 | 2888690 | 1.471e-3 | 0.305 | fail |

Strang vs Lie takes 0.722 s against 0.455 s (16 threads, `--repeat 3`). A/B
of the two other gate passes, alternating the ticket-03 build (`4a46fd0`) and
the exploration build, medians of per-run medians (16 threads: 4 runs of
`--repeat 5`; 1 thread: 3 runs of `--repeat 3`):

| threads | step-doubling [s] | average, K=1 [s] | S2, K=2/3 [s] |
|---|---|---|---|
| 16 | 0.428 | 0.430 (+0.4%) | 0.381 (-10.9%) |
| 1 | 2.515 | 2.317 (-7.9%) | 2.101 (-16.5%) |

The average gains nothing at 16 threads, probably because its 30% more steps
add barriers and serial work (diffusion, background integration). S2 with
K=2/3 passes the gate and the keep rule, but K was picked on this very case.
Two more end times of the same case were checked, against references written
the same way (`--t-end 1e-8 --write-reference` and `--t-end 1e-6
--write-reference`, 378 and 692 steps):

| t_end | estimator | K | steps | nfev | max total err | max scaled err |
|---|---|---|---|---|---|---|
| 1e-8 | step-doubling | - | 35 | 1719862 | 1.237e-3 | 0.232 |
| 1e-8 | reversed, S2 | 2/3 | 39 | 1320330 | 2.685e-3 | 0.533 |
| 1e-8 | reversed, S2 | 3 | 56 | 1887898 | 1.285e-3 | 0.272 |
| 1e-8 | reversed, average | 1 | 43 | 1486409 | 1.734e-3 | 0.280 |
| 1e-8 | reversed, average | 2 | 50 | 1730724 | 1.123e-3 | 0.186 |
| 1e-6 | step-doubling | - | 54 | 5727857 | 4.093e-4 | 0.0271 |
| 1e-6 | reversed, S2 | 2/3 | 63 | 4797499 | 1.021e-3 | 0.0500 |
| 1e-6 | reversed, S2 | 3 | 94 | 6411875 | 3.856e-4 | 0.0204 |
| 1e-6 | reversed, average | 1 | 70 | 5166518 | 1.471e-4 | 0.0198 |

At these end times, S2 with K=2/3 is 1.8 to 2.5 times less accurate than
step-doubling. To match step-doubling it needs K of about 3, which costs more
nfev than step-doubling. The average needs K of about 2 at 1e-8, again with no
gain. Dropped: the code is reverted, step-doubling and `strang_pair` stay, and
there is no ADR.

Why it loses: step-doubling keeps the two-half-step solution, a Strang
solution with step dt/2. Strang is second order, so that solution has about a
quarter of the error of one Strang step of dt. Rerunning step-doubling's own
steps but keeping a single step of dt shows this (max total / max scaled
error):

| t_end | two half-steps (kept) | S1 over dt | S2 over dt |
|---|---|---|---|
| 1e-8 | 1.237e-3 / 0.232 | 4.980e-3 / 0.943 | 4.056e-3 / 0.792 |
| 1e-7 | 1.496e-3 / 0.193 | 5.947e-3 / 0.784 | 1.986e-3 / 0.326 |
| 1e-6 | 4.093e-4 / 0.0271 | 1.627e-3 / 0.125 | 1.622e-3 / 0.0800 |

An embedded pair keeps one step of dt, so it needs about half the dt (twice
the steps) for the same accuracy. With 3 solves per attempt, both candidates
then do 6 solves where step-doubling does 5. Their estimates also fail in
different ways:

- Lie's error is O(dt^2) and much larger than Strang's here, so the steps get
  too small.
- The reversed pair is a Milne-type estimate of the wrong kind. To leading
  order, S1's error is `(2P - Q)/24 dt^3` and S2's is `(2Q - P)/24 dt^3`, with
  `P = [D,[D,R]]` and `Q = [R,[R,D]]`. Their difference `(P - Q)/8 dt^3` is
  not proportional to either error and cancels when P and Q are close. At
  1e-7, one S2 step happens to be nearly as accurate as the two half-steps
  (table above), which is why S2 passes the gate there but not at the other
  two end times.

Not tried: keep step-doubling and take the first half-reaction of both paths
from one solve (from `y0` to `t + dt/2` with an rsoda `tcrit` stop at
`t + dt/4`). That is 4 solves per attempt instead of 5 with the same
estimator, but results would no longer be bit-identical to now.
