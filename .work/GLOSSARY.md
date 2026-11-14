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

### Voxel model

**Voxel grid**:
The uniform Cartesian nx×ny×nz domain with per-axis spacing and zero-flux boundaries; each voxel holds one full species state vector.
_Avoid_: mesh, lattice, cell (cell means the Phase-1 unit cell)

**Seed**:
A set of Geant4-DNA molecules (CSV: species, x, y, z) binned into voxel concentrations and injected into the grid at a given time.
_Avoid_: track input, initial condition (a seed may arrive later than t0)

**Seed alias**:
The mapping from a Geant4-DNA molecule name (`GetName()` form like `HO_2°^0`, or the user tag like `HO2°`) to a radkin species label; built-in for every Geant4-DNA chemistry option3 solute, extendable by `[voxel.geant4.aliases]`.
_Avoid_: species map, translation

**Splitting step**:
One Strang step reaction(dt/2) → diffusion(dt) → reaction(dt/2) over the whole grid; its dt is adapted by step-doubling.
_Avoid_: time step (ambiguous with rsoda's internal steps)

**Snapshot**:
The full 3D concentration field of every species at a requested output time.
_Avoid_: frame, dump

**Background voxel**:
A voxel whose state still matches the shared background vector (the config's initial solution, no seed species) within a tolerance far below the solver atol; it is advanced by copying one shared background integration instead of its own rsoda run.
_Avoid_: empty voxel (it holds O2, H3O+, OH- in aerated runs), inactive voxel

**Track segment run**:
A Geant4-DNA run where only a small water region (2x2x10 um at the track entrance) holds the track, so its deposited energy, its molecules and its IRT G(t) all refer to that region; the fine-voxel radkin reference case.
_Avoid_: cropped track (cutting a cube from a whole-track run is a different, inexact thing)

### Performance

**Time breakdown**:
The split of one solve's wall time into kernel RHS, kernel Jacobian, linear-solver factor, linear-solver solve, and the remainder (rsoda bookkeeping and output storage), measured by the instrumented bench.
_Avoid_: profile (also means a sampling profile), timings

**Per-RHS cost**:
Solve time divided by the number of RHS evaluations (nfev), solver overhead included; the headline figure compared against radiopyo.
_Avoid_: RHS time (that is the kernel's share only)

**Coarsening level**:
The global resolution of the voxel grid relative to the configured fine grid: level L merges every 2^L×2^L×2^L block of fine voxels into one uniform voxel. The whole grid is always at one level.
_Avoid_: AMR level, refinement level, octree depth
