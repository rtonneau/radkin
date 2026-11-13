---
bump: minor
floor: minor
---
- Sparse Jacobian option: `--sparse-jacobian` (library: `RunOptions::jacobian = JacobianMode::Sparse`) solves the stiff steps with a sparse LU on the reaction network's sparsity pattern; results match the dense run within the golden tolerances.
- `examples/bench.rs --jacobian dense|sparse|both` times the dense and sparse solvers.
- rsoda updated to 0.2.0 (pluggable linear solvers, sparse LU behind its `sparse` feature).
