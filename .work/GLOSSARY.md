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

### Stochastic solver

**SSA**:
The exact stochastic simulation of the reaction network on integer molecule counts (Gillespie), selected by `solver = "ssa"`; the reference that tau-leaping and the LSODA comparison are checked against.
_Avoid_: Gillespie solver (also used loosely for tau-leaping), Monte Carlo

**Propensity**:
The probability per unit time that one reaction (or one hop) fires, computed from the molecule counts and the volume; the stochastic counterpart of a reaction rate.
_Avoid_: rate (means the mol/L/s rate of the LSODA path), hazard

**Trajectory**:
One stochastic realisation of a run, drawn from its own random stream (config seed, stream index = trajectory index).
_Avoid_: replica, sample, run (a run is the whole ensemble)

**Ensemble**:
The N trajectories of one run; outputs are its mean and standard deviation on the output time grid.
_Avoid_: batch, population

**Tau-leaping**:
The opt-in approximate stepping that fires many reactions (and hops) per step tau with Poisson counts, with critical reactions still fired one at a time (Cao-Gillespie-Petzold 2006).
_Avoid_: leaping SSA, approximate SSA

**Hop**:
The move of one molecule to a face-neighbouring voxel, an event of propensity `D/h²` per molecule and face; it replaces the diffusion step in stochastic voxel runs.
_Avoid_: jump, diffusion event

**Critical spacing**:
The voxel spacing `h* = 0.2527 κ / D` of a bimolecular reaction (`κ` its pair rate constant, `D` the summed diffusion coefficients) below which no propensity can reproduce its mean encounter time on the grid (Hellander 2012). For `2 A` reactions `κ = 2k`.
_Avoid_: validity bound, h_min (h_min is where the uncorrected RDME starts to degrade)

**Mesh correction**:
The scaling `1 / (1 - h*/h)` of a bimolecular propensity on a grid level of spacing `h > h*`, which makes the RDME's mean encounter time match the macroscopic rate constant (`[ssa] rate_correction = "hellander"`).
_Avoid_: rate correction (the config key only), mesoscopic rate, renormalisation

**Reversible pair**:
A forward/backward reaction pair, declared with `<->` in the network or inferred from two one-way mass-action reactions whose sides swap once fixed species are dropped (each the single candidate reverse of the other, neither a direction of a declared pair). In a voxel SSA run, its two directions share one mesh-correction factor, taken from the direction that is bimolecular once fixed species are dropped.
_Avoid_: equilibrium reaction (only one parsed kind of pair), reverse reaction

**Background run**:
The seedless RDME companion run, on the same grid and settings, whose counts are subtracted to give background-net yields.
_Avoid_: blank run, control run

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

**Handover time**:
The time at which the particle chemistry stops and its molecules are binned onto the voxel grid, which computes the chemistry from there on. The particle chemistry is either a Geant4-DNA run (IRT_syn, whose molecules become radkin seeds) or radkin's own IRT stage (one fixed time per run, in-process).
_Avoid_: start time (ambiguous with the config's t0), seed time (a seed time is any injection time), coupling time, switch time

### IRT stage

**IRT stage**:
The part of a run where every molecule is a particle with a position, and reactions are drawn as independent reaction times per pair (and per molecule for background reactions); it runs from the 1 ps molecules to the handover time, or to the end of the run.
_Avoid_: particle stage, spur stage, microscopic stage, IRT_syn (the name of Geant4's implementation only)

**Full resampling**:
The pair rule of the IRT stage that discards every sampled reaction time after each reaction and draws them all again from the current positions; the default and the reference rule.
_Avoid_: synchronous mode, IRT_syn mode (Geant4's IRT_syn keeps its reaction set between reaction steps)

**Incremental resampling**:
The pair rule of the IRT stage that keeps the sampled reaction times and draws new ones only for pairs involving the molecules that reacted, were created or were injected; the working rule for many tracks.
_Avoid_: classic IRT, lazy resampling

**Pair reaction**:
A reaction of the IRT stage between two particles, with a type (fully or partially diffusion-controlled), an effective radius and a sampled reaction time.
_Avoid_: bimolecular reaction (a reaction with a background species is bimolecular too, but is a molecule reaction here)

**Molecule reaction**:
A reaction of the IRT stage of one particle alone, with an exponential reaction time: first order, with a fixed species such as water, or with a background species.
_Avoid_: scavenging reaction (only one of its three cases), first-order reaction

**Background species**:
In the IRT stage, a species with an initial concentration (O2, H3O+ and OH- at the pH of the run): particles react with it at the rate `k c`, it keeps its concentration, and what it loses is counted, not debited. The molecules of the same species made by a track are particles.
_Avoid_: scavenger (O2 is one, H3O+ is not called one), bulk species

**Contact rule**:
A pair injected inside its effective radius reacts at once, with the contact probability of its reaction (1 when fully diffusion-controlled). It applies at injection only; Geant4 applies it at the start of its chemistry (1 ps).
_Avoid_: overlap reaction, zero-time reaction

**Image pair**:
A particle and the mirror image of its partner across a reflective wall, sampled as one more independent pair so that a pair near a wall reacts as often as in free space.
_Avoid_: ghost particle, periodic image (the walls reflect, the box is not periodic)

**Outflow box**:
The box of an IRT run where a molecule that a move puts outside is removed, as Geant4-DNA does at the faces of its water volume; used only to compare with a Geant4 run made in a small box.
_Avoid_: absorbing wall (it loses about half as many molecules), open boundary

### Performance

**Time breakdown**:
The split of one solve's wall time into kernel RHS, kernel Jacobian, linear-solver factor, linear-solver solve, and the remainder (rsoda bookkeeping and output storage), measured by the instrumented bench.
_Avoid_: profile (also means a sampling profile), timings

**Per-RHS cost**:
Solve time divided by the number of RHS evaluations (nfev), solver overhead included; the headline figure compared against radiopyo.
_Avoid_: RHS time (that is the kernel's share only)

**Coarsening level**:
The global resolution of the voxel grid relative to the configured fine grid: level L merges every 2^L×2^L×2^L block of fine voxels into one uniform voxel. The whole grid is at one coarsening level; under the RDME it is a floor for every species level.
_Avoid_: AMR level, refinement level, octree depth

**Species level**:
The fixed level of the 2^L hierarchy on which one species lives in a stochastic voxel run (`[ssa.species_levels]`, default 0); its molecules hop and react on that level's voxels, and a species sits at the larger of its species level and the current coarsening level.
_Avoid_: species grid, species resolution, multi-grid level
