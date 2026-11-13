---
bump: minor
floor: minor
---
- Solves are 16-25% faster on the benchmark cases, with identical results (same steps, RHS and Jacobian counts), from a faster right-hand-side kernel and rsoda 0.2.1's faster step loop.
- New `timing` cargo feature and `examples/bench.rs --breakdown`: the time a solve spends in the kernel RHS, kernel Jacobian, LU factorization, linear solves and the rest. `tools/breakdown.py` runs it on the benchmark matrix, and `docs/perf/` records the results.
