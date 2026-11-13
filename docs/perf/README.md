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

_Filled by the optimization tickets: each attempt, kept or dropped, with its
before/after plain medians._
