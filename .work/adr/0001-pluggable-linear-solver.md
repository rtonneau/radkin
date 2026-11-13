# rsoda linear solve behind a public trait, sparse LU as optional faer feature

_ADR 0001, recorded 2026-11-13 (session 2026-11-13__rsoda-sparse-jacobian)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

rsoda only had dense and banded LU, which cannot scale to Phase-2 voxel systems, and the Phase-2 choice between operator splitting and a monolithic solve is still open. We put the iteration-matrix assemble/factor/solve/norm behind a public linear-solver trait (dense and banded become bit-identical implementations of it) and ship a sparse LU built on faer behind an optional `sparse` cargo feature, so rsoda keeps zero dependencies by default and users can plug their own solvers (e.g. Krylov or block-diagonal) later.

Considered: built-in faer dependency (forces a dependency on every user), hand-written sparse LU (numerically delicate code to own), Krylov/LSODKR port (different convergence tests and preconditioning; much larger port, still possible through the trait).
