# radkin glossary

radkin is a Rust port of radiopyo: homogeneous chemical kinetics of water radiolysis (G-values, rate equations, beams), integrated with the LSODA solver rsoda, with a planned voxelised diffusion phase.

<!-- Format adapted from mattpocock-skills' domain-modeling (GLOSSARY-FORMAT.md): https://github.com/mattpocock/skills.git -->

<!--
How to add a term (the moment it is resolved, not in a batch at the end):

**Term**:
One or two sentences: what it IS, not what it does.
_Avoid_: the other words people use for it

- Be opinionated: one word per concept; list the others under _Avoid_.
- Only terms specific to this project. General programming concepts and
  implementation details don't belong here: this is a glossary, not a spec.
- Group terms under "### <area>" subheadings once natural clusters emerge.
-->

## Language

### Solver

**Iteration matrix**:
The Newton matrix `P = I - h*el0*J` that rsoda factors and solves in the stiff (BDF) corrector.
_Avoid_: Jacobian matrix (J is only one ingredient of P), Newton matrix

**Linear solver**:
The pluggable rsoda component that owns the iteration matrix's storage: it assembles P from J, factors it, solves with it, and reports the weighted norm of J for stiffness detection. Dense, banded and sparse are implementations.
_Avoid_: backend, matrix storage

**Sparsity pattern**:
The fixed set of (row, column) positions where J may be nonzero, given once before the run (CSC); every Jacobian evaluation writes values only at those positions.
_Avoid_: structure, mask

**Column coloring**:
A grouping of Jacobian columns that share no nonzero row, so one finite-difference RHS call per group estimates all of its columns (Curtis-Powell-Reid).
_Avoid_: column grouping, compression

### Performance

**Time breakdown**:
The split of one solve's wall time into kernel RHS, kernel Jacobian, linear-solver factor, linear-solver solve, and the remainder (rsoda bookkeeping and output storage), measured by the instrumented bench.
_Avoid_: profile (also means a sampling profile), timings

**Per-RHS cost**:
Solve time divided by the number of RHS evaluations (nfev), solver overhead included; the headline figure compared against radiopyo.
_Avoid_: RHS time (that is the kernel's share only)
