# radkin performance: time breakdown of a solve

This page records where a radkin solve spends its time and what was done about
it (gps session `2026-11-13__rsoda-per-step-overhead`, idea
`.work/idea/rsoda-per-step-overhead.md`).

## How it is measured

- `examples/bench.rs` solves a case `--repeat` times in-process (see
  `tools/benchmark.py`). Built with the `timing` cargo feature (which turns on
  rsoda's `timing` feature), `--breakdown` adds, per repeat, the time spent in:
  - **rhs**: kernel right-hand-side calls (`kernel::derive`);
  - **jac**: kernel Jacobian calls (dense or sparse);
  - **factor**: rsoda `LinearSolver::prepare` minus **jac**, i.e. assembly and
    LU factorization of the iteration matrix `P = I - h*el0*J`;
  - **lin**: rsoda `LinearSolver::solve` calls (solves with the factors);
  - **rest**: the instrumented solve time minus the four above: rsoda's step
    bookkeeping (prediction, corrector loop, error norms, step and order
    selection, method switching) and radkin's per-step output storage and
    per-chunk setup.
- `tools/breakdown.py` runs every case of the benchmark matrix, dense and
  sparse, with a plain build (median solve time, no timers) and a `timing`
  build (`target/timing`, the shares). Shares are relative to the instrumented
  median; **timer** is the instrumented median over the plain median.
- **Per-RHS cost** is the plain median solve time divided by nfev.
- Keep rule for optimizations: a change stays only if it cuts the plain
  median solve time by at least 5% on the cases it targets, with no case
  slower beyond noise, and every golden test passes at unchanged tolerances.

Reproduce:

```bash
python -I tools/breakdown.py --repeat 10 --out breakdown.json
```

## Baseline

radkin `52c0e47`, rsoda `08a2abf` (0.2.0 + `timing`). Intel Core i7-7820X
3.6 GHz, Windows 11, rustc 1.99.1, `--repeat 10`, medians. n_dyn = 22 for
Labarbe.

| case | jac | solve [ms] | nfev | nst | njev | per-RHS [us] | rhs | jac | factor | lin | rest | timer |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| labarbe run 1e-7..10 s | dense | 7.494 | 4298 | 2567 | 158 | 1.74 | 35% | 4% | 7% | 15% | 38% | x1.04 |
| labarbe run 1e-7..10 s | sparse | 8.088 | 4316 | 2573 | 138 | 1.87 | 33% | 4% | 11% | 14% | 38% | x1.06 |
| labarbe run, 200 t_eval | dense | 7.065 | 4298 | 2567 | 158 | 1.64 | 37% | 4% | 7% | 15% | 37% | x1.09 |
| labarbe run, 200 t_eval | sparse | 7.674 | 4316 | 2573 | 138 | 1.78 | 34% | 4% | 12% | 15% | 35% | x1.06 |
| water constant run 0..1 s | dense | 4.050 | 2568 | 1771 | 200 | 1.58 | 30% | 8% | 8% | 15% | 39% | x1.07 |
| water constant run 0..1 s | sparse | 5.249 | 2706 | 1810 | 217 | 1.94 | 24% | 7% | 20% | 16% | 32% | x1.13 |
| water single pulse chunked | dense | 4.242 | 3085 | 2028 | 120 | 1.38 | 34% | 5% | 5% | 13% | 44% | x1.08 |
| water single pulse chunked | sparse | 5.573 | 3293 | 2140 | 150 | 1.69 | 29% | 5% | 14% | 15% | 37% | x1.07 |
| water pulsed 100 Hz chunked 0..1 s | dense | 148.707 | 122551 | 75379 | 3899 | 1.21 | 39% | 4% | 4% | 13% | 41% | x1.08 |
| water pulsed 100 Hz chunked 0..1 s | sparse | 178.302 | 122421 | 75285 | 3853 | 1.46 | 32% | 4% | 11% | 15% | 38% | x1.09 |
| water macro-pulsed chunked 0..0.1 s | dense | 8.207 | 6519 | 4443 | 224 | 1.26 | 37% | 4% | 4% | 12% | 43% | x1.11 |
| water macro-pulsed chunked 0..0.1 s | sparse | 10.678 | 6647 | 4489 | 233 | 1.61 | 29% | 4% | 11% | 13% | 43% | x1.09 |
| labarbe pulsed run 0..0.015 s | dense | 37.577 | 18815 | 17041 | 916 | 2.00 | 31% | 5% | 8% | 16% | 40% | x1.05 |
| labarbe pulsed run 0..0.015 s | sparse | 41.050 | 18977 | 17102 | 904 | 2.16 | 28% | 5% | 14% | 17% | 37% | x1.07 |
| paper single pulse 0..1 s | dense | 4.875 | 4078 | 2326 | 113 | 1.20 | 40% | 4% | 4% | 9% | 43% | x1.02 |
| paper single pulse 0..1 s | sparse | 5.677 | 4096 | 2331 | 109 | 1.39 | 35% | 4% | 10% | 11% | 40% | x1.04 |
| paper pulsed chunked 0..0.1 s | dense | 14.421 | 12158 | 7126 | 363 | 1.19 | 39% | 4% | 4% | 11% | 42% | x1.09 |
| paper pulsed chunked 0..0.1 s | sparse | 17.429 | 12642 | 7271 | 343 | 1.38 | 35% | 4% | 10% | 13% | 38% | x1.10 |

Reading:

- **rest** is the largest share on every case (32-44%). This is the rsoda
  step loop and radkin's output, not the linear algebra.
- **rhs** (kernel derive, about 0.5 us per call) is the next share (24-40%).
- **lin** (dense 22x22 forward/back substitution, about 0.4 us per call,
  called on every corrector iteration) is 9-17%. **factor** is small for
  dense (4-8%) and larger for sparse (10-20%).
- **jac** is small (4-8%): the Jacobian is reused for many steps (njev is
  about nst/20).
- The timers cost 2-13% of the solve (Instant::now around each sub-µs call),
  so the shares are approximate. Decisions use the plain medians.

## Optimizations

Each attempt, kept or dropped, with its before/after plain medians. A/B
numbers come from the two bench binaries run in alternation on the same
machine (6 rounds x 10 repeats per case, medians).

### Kernel right-hand side (ticket 03), kept

`derive` spent about 250 ns in the rates and 250 ns in the species balance
(Labarbe: 85 reactions, 1-4 reactant slots, 1-5 stoichiometry entries).
Microbenchmarks of `derive` alone (Labarbe state, best of 3):

| variant | derive [ns] | result |
|---|---|---|
| baseline | ~520 | |
| skip `sw`/`sl` when unused | ~515 | -1%, kept as part of the final change |
| flat CSR arrays for slots and stoichiometry, fused rate + scatter | ~547 | slower, dropped |
| species-major gather (transposed stoichiometry) | ~600 | slower (long add chains on frequently hit species), dropped |
| fixed-width packed reactions (3 slots, 4 rows, ordered tail) | ~485 | kept |
| packed + unchecked indexing | ~470 | bounds checks are only ~3%; `unsafe` not worth it, dropped |
| packed + skip `sw`/`sl` (final) | ~440 | -15% |

The final change, bit-identical to the old `derive` (unit tests compare the
bits, also with synthetic negative products): reactions are packed with
fixed-width reactant slots (padded with the factor-1 species) and fixed-width
stoichiometry rows (padded with a dummy accumulator row of coefficient 0), so
the inner loops have a fixed trip count. Longer rows continue in a tail, in
order. The depletion switches `sw` are computed only when a reaction has a
negative product or a species is below 0, and their slopes `sl` only for
Jacobians. The reordered sums (gather, split accumulators) were not pursued:
the balances cancel heavily, and the RHS golden check is at 1e-12.

Solve times (plain build), nfev and nst unchanged on every case:

| case | jac | old [ms] | new [ms] | change | nfev old/new | nst old/new |
|---|---|---|---|---|---|---|
| labarbe run 1e-7..10 s | dense | 7.679 | 7.106 | -7.5% | 4298/4298 | 2567/2567 |
| labarbe run 1e-7..10 s | sparse | 8.126 | 7.529 | -7.3% | 4316/4316 | 2573/2573 |
| labarbe run, 200 t_eval | dense | 7.211 | 6.677 | -7.4% | 4298/4298 | 2567/2567 |
| labarbe run, 200 t_eval | sparse | 7.687 | 7.068 | -8.1% | 4316/4316 | 2573/2573 |
| water constant run 0..1 s | dense | 4.176 | 3.788 | -9.3% | 2568/2568 | 1771/1771 |
| water constant run 0..1 s | sparse | 5.567 | 5.110 | -8.2% | 2706/2706 | 1810/1810 |
| water single pulse chunked | dense | 4.214 | 3.915 | -7.1% | 3085/3085 | 2028/2028 |
| water single pulse chunked | sparse | 5.561 | 5.174 | -7.0% | 3293/3293 | 2140/2140 |
| water pulsed 100 Hz chunked 0..1 s | dense | 148.053 | 133.663 | -9.7% | 122551/122551 | 75379/75379 |
| water pulsed 100 Hz chunked 0..1 s | sparse | 176.795 | 165.301 | -6.5% | 122421/122421 | 75285/75285 |
| water macro-pulsed chunked 0..0.1 s | dense | 8.158 | 7.621 | -6.6% | 6519/6519 | 4443/4443 |
| water macro-pulsed chunked 0..0.1 s | sparse | 10.767 | 9.889 | -8.2% | 6647/6647 | 4489/4489 |
| labarbe pulsed run 0..0.015 s | dense | 37.477 | 35.940 | -4.1% | 18815/18815 | 17041/17041 |
| labarbe pulsed run 0..0.015 s | sparse | 41.148 | 38.759 | -5.8% | 18977/18977 | 17102/17102 |
| paper single pulse 0..1 s | dense | 4.941 | 4.531 | -8.3% | 4078/4078 | 2326/2326 |
| paper single pulse 0..1 s | sparse | 5.720 | 5.240 | -8.4% | 4096/4096 | 2331/2331 |
| paper pulsed chunked 0..0.1 s | dense | 14.559 | 13.246 | -9.0% | 12158/12158 | 7126/7126 |
| paper pulsed chunked 0..0.1 s | sparse | 17.581 | 16.015 | -8.9% | 12642/12642 | 7271/7271 |

All cases are faster (-4.1% to -9.7%, median about -8%).

The dense and sparse kernel Jacobians are 4-8% of the solve. They were left
as they are.

### Dense LU factor and solve in rsoda (ticket 04), dropped

The dense `dgefa`/`dgesl` inner loops were rewritten on column slices (the
same element operations in the same order, bit-identical, no bounds checks).
Microbenchmark, 22x22 matrix of 30% density: `dgefa` 2.9 -> 1.9 us (-33%),
`dgesl` 360 -> 300 ns (-15%). `dgesl` is close to its latency bound: back
substitution chains one division per row, and removing it (multiplying by
stored reciprocals) would change the results.

At the solve level the gain is below the 5% rule. Dense cases gain -1.2% to
-4.8%, around -3% (factor + lin are about 20% of a dense solve). The sparse
cases, which do not use this code, moved by +1% to +8%, which gives the
noise floor between two builds:

{ab}

Not kept. The patch is in the ticket-04 log of the session
(`.work/sessions/2026-11-13__rsoda-per-step-overhead/03-implement/04-rsoda-dense-lu/`)
in case a later need (e.g. larger dense systems) makes it pay.

### rsoda step bookkeeping (ticket 05), kept

A temporary cycle-counter probe (`rdtsc` per DSTODA jump target, not
committed) split rsoda's run on Labarbe (dense) as: RHS calls 32%, linear
solves 15%, prepare 13%, and rsoda's own work about 32%. The latter was
predict (9%), the corrector's vector updates (about 9%), the history update
after a step (7%), and the driver's error weights and checks (3%). These were
indexed loops over the history array and the work vectors, with bounds checks
and no vectorization.

Kept: `predict`/`retract` run as column updates on slices (each sweep is
`col(c) += col(c + 1)` for ascending columns, the same values as the indexed
sweep, checked bit for bit by a unit test), and the corrector, history update
and rescale loops zip their vectors. Same element operations in the same
order, so the results are bit-identical. rsoda's golden counters and states
are unchanged, and nfev/nst are unchanged on every case below.

{ab}

All cases are faster (-5.8% to -18.4%, dense around -15%, sparse around
-11%).

Dropped:

- **4-lane `vmnorm`** (four running maxima, bit-identical because the max
  ignores NaN and does not depend on order): -0.6% to -3.6%, below the rule.

{vm}

- **Flat output storage in `solver::solve`** (one `Vec<f64>` instead of a
  `Vec` per stored step). The upper bound, not storing at all, is -6% to
  -8% on run-mode cases. A flat buffer that still copies each step gives 0%
  to -7.6% (most cases -1% to -3.6%). It would also change the public
  `OdeResult::y()` layout. Not worth it.

## Final

radkin pinned to rsoda `f5d1312` (0.2.1: `timing` feature and the step-loop
slices, rsoda PR #4). Baseline (radkin `6dc10f5`, rsoda `08a2abf`) against
final, same A/B protocol (6 rounds x 10 repeats, plain builds):

| case | jac | old [ms] | new [ms] | change | nfev old/new | nst old/new |
|---|---|---|---|---|---|---|
| labarbe run 1e-7..10 s | dense | 7.849 | 5.930 | -24.4% | 4298/4298 | 2567/2567 |
| labarbe run 1e-7..10 s | sparse | 8.378 | 6.337 | -24.4% | 4316/4316 | 2573/2573 |
| labarbe run, 200 t_eval | dense | 7.202 | 5.374 | -25.4% | 4298/4298 | 2567/2567 |
| labarbe run, 200 t_eval | sparse | 7.668 | 5.747 | -25.0% | 4316/4316 | 2573/2573 |
| water constant run 0..1 s | dense | 3.957 | 3.261 | -17.6% | 2568/2568 | 1771/1771 |
| water constant run 0..1 s | sparse | 5.335 | 4.449 | -16.6% | 2706/2706 | 1810/1810 |
| water single pulse chunked | dense | 4.161 | 3.241 | -22.1% | 3085/3085 | 2028/2028 |
| water single pulse chunked | sparse | 5.437 | 4.538 | -16.5% | 3293/3293 | 2140/2140 |
| water pulsed 100 Hz chunked 0..1 s | dense | 143.722 | 112.992 | -21.4% | 122551/122551 | 75379/75379 |
| water pulsed 100 Hz chunked 0..1 s | sparse | 175.615 | 141.626 | -19.4% | 122421/122421 | 75285/75285 |
| water macro-pulsed chunked 0..0.1 s | dense | 8.206 | 6.436 | -21.6% | 6519/6519 | 4443/4443 |
| water macro-pulsed chunked 0..0.1 s | sparse | 10.829 | 8.690 | -19.8% | 6647/6647 | 4489/4489 |
| labarbe pulsed run 0..0.015 s | dense | 37.325 | 28.595 | -23.4% | 18815/18815 | 17041/17041 |
| labarbe pulsed run 0..0.015 s | sparse | 40.579 | 31.650 | -22.0% | 18977/18977 | 17102/17102 |
| paper single pulse 0..1 s | dense | 4.903 | 3.863 | -21.2% | 4078/4078 | 2326/2326 |
| paper single pulse 0..1 s | sparse | 5.696 | 4.541 | -20.3% | 4096/4096 | 2331/2331 |
| paper pulsed chunked 0..0.1 s | dense | 14.318 | 11.334 | -20.8% | 12158/12158 | 7126/7126 |
| paper pulsed chunked 0..0.1 s | sparse | 17.521 | 14.022 | -20.0% | 12642/12642 | 7271/7271 |

Every case is 16.5% to 25.4% faster (dense about -22%, sparse about -20%),
with the same nfev and nst, so the step sequence is unchanged. All changes
kept are bit-identical to the code they replace.

Final breakdown (`tools/breakdown.py --repeat 10`):

| case | jac | solve [ms] | nfev | nst | njev | per-RHS [us] | rhs | jac | factor | lin | rest | timer |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| labarbe run 1e-7..10 s | dense | 5.808 | 4298 | 2567 | 158 | 1.35 | 36% | 6% | 9% | 18% | 31% | x1.04 |
| labarbe run 1e-7..10 s | sparse | 6.193 | 4316 | 2573 | 138 | 1.43 | 33% | 5% | 14% | 18% | 30% | x1.10 |
| labarbe run, 200 t_eval | dense | 5.617 | 4298 | 2567 | 158 | 1.31 | 38% | 6% | 9% | 19% | 28% | x1.05 |
| labarbe run, 200 t_eval | sparse | 5.794 | 4316 | 2573 | 138 | 1.34 | 35% | 5% | 15% | 19% | 27% | x1.10 |
| water constant run 0..1 s | dense | 3.301 | 2568 | 1771 | 200 | 1.29 | 30% | 10% | 9% | 18% | 33% | x1.07 |
| water constant run 0..1 s | sparse | 4.574 | 2706 | 1810 | 217 | 1.69 | 23% | 9% | 23% | 19% | 26% | x1.07 |
| water single pulse chunked | dense | 3.458 | 3085 | 2028 | 120 | 1.12 | 34% | 6% | 5% | 17% | 38% | x1.04 |
| water single pulse chunked | sparse | 4.689 | 3293 | 2140 | 150 | 1.42 | 27% | 6% | 16% | 18% | 33% | x1.08 |
| water pulsed 100 Hz chunked 0..1 s | dense | 118.299 | 122551 | 75379 | 3899 | 0.97 | 39% | 5% | 5% | 16% | 35% | x1.07 |
| water pulsed 100 Hz chunked 0..1 s | sparse | 148.499 | 122421 | 75285 | 3853 | 1.21 | 32% | 5% | 13% | 18% | 33% | x1.08 |
| water macro-pulsed chunked 0..0.1 s | dense | 6.764 | 6519 | 4443 | 224 | 1.04 | 37% | 6% | 5% | 14% | 38% | x1.07 |
| water macro-pulsed chunked 0..0.1 s | sparse | 9.061 | 6647 | 4489 | 233 | 1.36 | 28% | 5% | 13% | 15% | 39% | x1.09 |
| labarbe pulsed run 0..0.015 s | dense | 29.976 | 18815 | 17041 | 916 | 1.59 | 31% | 6% | 10% | 20% | 33% | x1.04 |
| labarbe pulsed run 0..0.015 s | sparse | 33.901 | 18977 | 17102 | 904 | 1.79 | 28% | 6% | 16% | 20% | 29% | x1.05 |
| paper single pulse 0..1 s | dense | 3.873 | 4078 | 2326 | 113 | 0.95 | 40% | 5% | 4% | 11% | 40% | x1.10 |
| paper single pulse 0..1 s | sparse | 4.643 | 4096 | 2331 | 109 | 1.13 | 35% | 5% | 12% | 13% | 35% | x1.06 |
| paper pulsed chunked 0..0.1 s | dense | 11.734 | 12158 | 7126 | 363 | 0.97 | 39% | 5% | 5% | 14% | 37% | x1.12 |
| paper pulsed chunked 0..0.1 s | sparse | 14.518 | 12642 | 7271 | 343 | 1.15 | 33% | 5% | 12% | 16% | 34% | x1.07 |

What is left: the kernel RHS (23-40%) and rsoda's step loop (rest, 26-40%)
are still the two largest shares, followed by the linear solves (11-20%).
None of the remaining candidates passed the 5% rule: dense LU slice loops
(about -3%), a 4-lane `vmnorm` (about -2.5%), flat output storage (about
-2%). Faster solves would need changes that are not bit-identical
(reordered sums in the kernel, reciprocal pivots in `dgesl`) or a different
algorithm (e.g. the Phase-2 per-voxel structure).
