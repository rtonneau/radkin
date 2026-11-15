# Stochastic solver (SSA)

`solver = "ssa"` simulates the reaction network on integer molecule counts with
Gillespie's stochastic simulation algorithm. The code is in `src/ssa/`; usage is in the
README. This file records how the solver is built and how it was validated.

## Model

The propensity model (`src/ssa/model.rs`) is built from the same `SimMatrices` as the
LSODA kernel, for a volume `V`. With `N_A V` the number of molecules per mol/L:

| Reaction | Propensity [1/s] |
|---|---|
| Mass action, rate `k · Π[reactants]` | `k · Π[fixed] · N_A V · Π n(n−1)…(n−m+1) / (N_A V)^m` over the dynamic reactants |
| `2 A -> …` (rate `k [A]²`, stoichiometry −2) | `k / (N_A V) · n (n − 1)`, two molecules removed per firing |
| Michaelis–Menten | its rate law at the concentrations `n / (N_A V)`, times `N_A V` |
| Radiolytic production (G-value) | `yield · dose rate · N_A V`, one molecule produced |

- `k` is the rate constant of the rate law, as stored in `SimMatrices`. For a scalar
  `k_value` of a single-reactant reaction (`2 A -> …`) radiopyo's convention divides the
  configured value by the stoichiometry when the config is read, so the propensity model
  takes the stored value as it is.
- Fixed species (H2O, or a buffered pH) enter as constant concentrations.
- A reaction fires only if every species it consumes on balance is present in sufficient
  number. This replaces the smooth depletion switch `c/(c+ε)` that the LSODA kernel puts
  on "negative" products.

## Exact SSA

`src/ssa/exact.rs` is Gillespie's direct method. The dose rate is piecewise constant
between the beam transitions (`src/ssa/schedule.rs`). When a segment or an output
interval ends, the pending event is discarded and redrawn with the new propensities;
waiting times are memoryless, so this is exact.

## Random numbers

One ChaCha8 generator per trajectory (`src/ssa/rng.rs`): the seed keys the generator and
the trajectory index selects its stream. The samplers use IEEE basic operations only
(the logarithm is `det_ln`, not the platform's libm), so a seed gives the same
trajectories on every platform. Ensemble statistics are exact integer sums accumulated
in trajectory order, so the output does not depend on the number of threads.

## Validation against LSODA

`tests/ssa_validation.rs` compares the ensemble mean with LSODA on the pure-water
chemistry (`data/advanced_water_radiolysis.toml`, 15 dynamic species), at 12 times
taken among the LSODA points after 1 µs:

| Case | Setup |
|---|---|
| Reaction-only | radicals at µmol/L level as initial state, no beam, 1 ms |
| Continuous beam | 1e4 Gy/s for 1 ms (10 Gy) |
| Pulsed beam | 3 pulses of 3 Gy (2 µs each) every millisecond, 3.5 ms |

Criterion, per species and time, in molecules:

```text
|mean − lsoda| ≤ 4 · SEM_eff + 1
```

`SEM_eff` is the larger of the sample standard error of the mean and the Poisson one,
`sqrt(lsoda / N)`. The second matters for species with a mean far below one molecule,
whose sample SEM is often 0. The seed is fixed, so the test is deterministic.

Measured with `V = 1e-14 L` (6022 molecules per µmol/L), 64 trajectories, seed 2026:

| Case | Checks | Worst deviation | Events per trajectory |
|---|---|---|---|
| Reaction-only | 180 | 2.0 SEM | 79 734 |
| Continuous beam | 180 | 2.8 SEM | 105 583 |
| Pulsed beam | 180 | 2.6 SEM | 164 860 |

The three cases run in 0.5 s in total. A factor 2 on the `2 A` propensity was tried by
hand: all three cases then fail, H2O2 being 30 SEM away from LSODA.

### Finite-volume difference

The ensemble mean of the stochastic model is not exactly the solution of the rate
equations. With 2000 to 4000 trajectories the only visible difference is on H3O+ and
OH−, in the two beam cases: their mean is about 4 molecules above LSODA.

| Volume | Trajectories | OH− LSODA | OH− SSA mean | Difference | SEM |
|---|---|---|---|---|---|
| 1e-15 L | 4000 | 59.8 | 64.3 | +4.5 molecules | 0.39 |
| 1e-14 L | 2000 | 598.4 | 602.6 | +4.2 molecules | 1.7 |

(pulsed case, t = 3.5 ms.) The difference is the same number of molecules at both
volumes, so its relative size falls as 1/V. It is expected: the G-values of H3O+ and of
e_aq are independent sources, so the net charge of a trajectory fluctuates, and the
water equilibrium `[H3O+][OH−] = Kw` is a convex function of that charge. The validation
volume is chosen so that this difference stays below half a SEM. In the reaction-only
case no such difference is seen (largest deviation 3.0 SEM over 180 checks at 4000
trajectories).

## Tau-leaping

`tau_leap = true` (`src/ssa/tau.rs`) replaces the exact steps by leaps in which every
reaction fires a Poisson number of times. It is the explicit tau-leaping of Cao,
Gillespie and Petzold:

- The leap bounds the relative change of every propensity by `tau_epsilon` (default
  0.03), from the mean and variance of the change of each reactant species (J. Chem.
  Phys. 124, 044109, 2006).
- A reaction with fewer than 10 firings left before it exhausts a species is critical:
  it fires one at a time, as in the exact SSA (J. Chem. Phys. 123, 054104, 2005).
- When a leap would fire fewer than 10 events, 100 exact steps are taken instead.
- A leap that would make a count negative is rejected and redone with half the step.
- A leap never crosses a beam transition or an output time.

Tau-leaping is an approximation. The reactions of a leap fire with the propensities of
its start, so the ensemble mean has a systematic error proportional to `tau_epsilon`, on
top of the sampling noise. Measured on the `radicals` case of
`examples/ssa_bench.rs` (three reactions, 1e-12 L, 600 000 to 1 800 000 molecules, 32
trajectories), as the largest difference of the two ensemble means:

| `tau_epsilon` | Difference [combined SEM] | Speed-up over exact |
|---|---|---|
| 0.03 | 32 | 330x to 650x |
| 0.01 | 11 | 220x |
| 0.003 | 3.7 | 80x |
| 0.001 | 2.4 | 25x |

At `tau_epsilon = 0.03` the difference is 0.36 % of the peak concentration of the
species. `tests/ssa_tau.rs` checks both parts with a fixed seed (1e-13 L, 64
trajectories): with `tau_epsilon = 0.001` the two ensembles agree within 4 combined SEM
(measured 2.4), and with 0.03 within 4 combined SEM plus 3 % of the species peak. The
standard deviations agree within their sampling error in both.

Speed-up on one core (`RAYON_NUM_THREADS=1`, `tau_epsilon = 0.03`):

| Case | Volume | Events per trajectory | Exact | Tau-leap | Speed-up |
|---|---|---|---|---|---|
| `radicals` | 1e-12 L | 2.0e6 | 0.142 s | 0.0002 s | 650x |
| `radicals` | 1e-13 L | 2.0e5 | 0.0138 s | 0.00025 s | 55x |
| `radicals` | 1e-14 L | 2.0e4 | 0.0014 s | 0.0003 s | 4.5x |
| `water-const` | 1e-13 L | 1.1e6 | 0.22 s | 0.25 s | 0.9x |
| `water-pulsed` | 1e-13 L | 1.6e6 | 0.32 s | 0.28 s | 1.1x |

(times per trajectory.) The speed-up grows with the molecule counts when no fast
equilibrium is present. On the water chemistry it is absent: the acid-base equilibria
keep a few species at very low counts with a very high turnover (at 1e-13 L, 2
molecules of HO2− react 1.4e8 times per second with water). The reactions that consume
them are critical and fire one at a time, which leaves about 17 events per leap. This
limit of explicit tau-leaping does not go away at larger volumes (1.1x at 1e-11 L).

## Voxel grid: the RDME

`radkin voxel` with `solver = "ssa"` solves the reaction-diffusion master equation
(`src/ssa/rdme.rs`, ADR 0003). Every voxel holds integer counts. Reactions fire inside
a voxel with the propensities above for the voxel volume. Diffusion is a hop of one
molecule to a face neighbour, with propensity `D / h²` per molecule and face;
boundary voxels have fewer faces (zero flux). The mean of the hop process obeys the same
discrete diffusion operator as the deterministic path.

- **Algorithm:** next-subvolume method. Each voxel keeps the time of its next event in
  an indexed binary heap; the earliest voxel fires one reaction or hop, and only the
  voxels whose counts changed draw a new time. A hop updates the reaction rate of its
  two voxels incrementally, from the reactions that read the hopping species.
- **Initial solution:** a Poisson draw per voxel and species with mean `c N_A V`.
- **Seeds:** the Geant4 molecules are injected as the integer counts they were binned
  from. Voxel runs have no beam, as on the deterministic path.
- **Ensemble:** `[ssa] trajectories` trajectories, trajectory `i` on stream `i` of
  `seed`, reduced in trajectory order. The outputs keep the deterministic layout with the
  ensemble mean, plus the standard deviation: `_std` snapshot files and `N_<sp>_std`,
  `G_<sp>_std` columns in `yields.csv`. The totals are recorded at 100 log-spaced times
  and at the snapshot and seed times.
- **Not available:** checkpoints and `--resume`.

### Validation against the deterministic voxel run

`tests/ssa_rdme.rs`, fixed seed, criterion `|mean − lsoda| ≤ 4 SEM_eff` plus the
tolerance of the deterministic run:

| Case | Compared | Worst deviation |
|---|---|---|
| Diffusion only: 4000 e_aq seeded in one corner voxel of an 8x2x2 grid of 50 nm voxels, 64 trajectories | mean field of every voxel at 3 times | 2.4 SEM |
| Pure water, 12 000 radicals seeded in the centre voxel of a 3x3x3 grid of 500 nm voxels, 48 trajectories | domain totals of every species at 4 times | 1.5 SEM |

The second case is a large-count grid (75 molecules per µmol/L and voxel): there the
RDME and the Strang-split rate equations agree.

### Coarsening

`voxel.max_coarsen_levels > 0` lets the RDME grid drop levels, as on the deterministic
path: the whole grid merges every 2x2x2 block into one voxel
(`src/ssa/rdme_coarsen.rs`).

- **Test.** The deterministic test (spread of the 8 children below a tolerance) can
  never pass on counts, because Poisson noise alone exceeds it. The grid drops a level
  when the counts of every block are consistent with a uniform distribution over its 8
  voxels. For a block holding `S` molecules of a species, Pearson's statistic
  `X = Σ (c_i − S/8)² / (S/8)` has the exact mean 7 and variance `14 (S − 1) / S` under a
  uniform split. For every species, two conditions must hold:
  - over all blocks with `S ≥ 1`: `(Σ X − 7 B) / sqrt(Σ variance) ≤ z`;
  - for each block with `S ≥ 16`: `(X − 7) / sqrt(variance) ≤ z + sqrt(2 ln B16)`, so
    that a localized cluster is not diluted by a large uniform background.
- **Threshold.** `z = 30 · coarsen_tol_factor`: 3 at the default factor 0.1.
- **Cadence.** The test runs every `8 × (number of voxels)` events.
- **Merge.** The 8 counts are summed. Hop rates and propensities use the spacing and
  volume of the current level.
- **Refine.** Before a seed is injected the grid returns to level 0: each molecule of a
  voxel goes to one of its 8 children at random. It may coarsen again afterwards.
- **Outputs** stay on the fine grid: a voxel of level `L` spreads its count evenly over
  its `8^L` fine voxels.

Species totals are conserved exactly through merging and refining.

Measured on a homogenising case (pure water, 4500 radicals seeded in the centre of an
8x8x8 grid of 100 nm voxels, 100 µs, 3 levels, 96 trajectories), against the fine-grid
ensemble, as the largest difference of the species totals over all record times:

| `coarsen_tol_factor` | `z` | Hops per trajectory | Fewer hops | Wall time | Difference [combined SEM] |
|---|---|---|---|---|---|
| fine grid | | 202 499 | | 0.97 s | |
| 0.03 | 0.9 | 52 998 | 3.8x | 0.35 s | 1.6 |
| 0.1 (default) | 3 | 28 534 | 7.1x | 0.20 s | 1.3 |
| 0.3 | 9 | 19 085 | 10.6x | 0.15 s | 2.6 |
| 1 | 30 | 12 528 | 16x | 0.12 s | 6.3 |
| 3 | 90 | 6 924 | 29x | 0.09 s | 24 |

Up to a factor of 0.3 the coarsened ensemble cannot be told from the fine one; from 1
on the grid coarsens while the radicals are still clustered and their recombination is
underestimated (more OH survives). The default 0.1 is kept. `tests/ssa_rdme.rs` checks
the default against the fine grid with a fixed seed (two seeds, so that the refinement is
exercised: 1.4 combined SEM, 5 times fewer hops).

### Tau-leaping on the grid

`tau_leap = true` also applies to voxel runs (`src/ssa/rdme_tau.rs`). The channels of a
leap are the reactions of every voxel and the hops of every species through every face
(`X(v) -> X(w)`, a first-order channel of propensity `n D / h²`).

- A channel is critical when the population it consumes is below 10. Critical channels
  fire one at a time.
- The leap is the Cao–Gillespie–Petzold bound over every (voxel, species) population,
  with the hops in and out counted in its mean and variance.
- Non-critical channels fire a Poisson number of times. A leap that would make a count
  negative is rejected and halved, so counts stay non-negative and hops conserve the
  molecules.
- A leap reads every voxel, so it is only taken when it replaces more than
  `10 + 4 × voxels` non-critical events. Otherwise the run takes exact next-subvolume
  events, in batches that double (from 8 to 512 events per voxel) while leaping keeps
  being refused.
- A leap never crosses a seed or snapshot time. The coarsening test counts the leaped
  events.

Measured against the exact RDME (three-reaction chemistry without fast equilibrium,
4x4x4 grid, 1 ms, 16 trajectories, `tau_epsilon = 0.03`):

| Case | Events per trajectory | Exact | Tau-leap | Speed-up | Difference of the totals |
|---|---|---|---|---|---|
| 1 µm voxels, 3000 OH per voxel, 60 000 e_aq and 60 000 H seeded in two corners | 4.1e6 | 2.08 s | 1.28 s | 1.6x | 0.12 % of the species peak |
| 1 µm voxels, 3000 molecules per voxel of three species, uniform | 1.0e7 | 5.70 s | 1.12 s | 5.1x | 0.31 % |
| 3 µm voxels, 81 000 molecules per voxel of three species, uniform | 3.9e7 | 19.96 s | 0.53 s | 38x | 0.25 % |

(wall time of the ensemble on 16 threads.) The first case is the fixed-seed test of
`tests/ssa_rdme.rs`: the totals agree within 4 combined SEM plus `tau_epsilon` times the
species peak, and the mean field of every voxel likewise. Leaping pays on smooth,
well-populated grids. A sharp front limits it: an empty voxel next to a full one may
only gain about one molecule per leap. On a radiolysis track at nm voxels every channel
is critical, no leap is taken, and the output is identical to the exact run; the leap
decisions then cost about 20 % more time (`data/voxel_example.toml`: 1.0 s exact, 1.2 s
with `tau_leap = true`, 2000 trajectories).

### Small voxels

At small voxels the plain RDME is biased. Two molecules can only react once they share
a voxel, and on a fine grid finding the same voxel takes longer than the macroscopic
encounter time. Each bimolecular reaction has a **critical spacing**

```text
h* = C₃ κ / (6 D) ≈ 0.2527 κ / D
```

where `κ` [m³/s] is the pair rate constant the RDME sees and `D` the summed
diffusion coefficient of the pair. `C₃ = 1.5164` is the mean number of lattice steps
to reach a target, divided by the number of voxels (Hellander, Hellander and Petzold,
Phys. Rev. E 85, 042901, 2012). `κ = k / N_A` for `A + B`; `κ = 2k / N_A` for `2 A`,
whose propensity `k/(N_A V) n (n − 1)` is `2k/(N_A V)` per pair, with `D = 2 D_A`.
radiopyo stores `k = k_value / 2` for `2 A`, so there `κ = k_value / N_A`. A fixed
reactant enters through `κ`: `e_aq + O_minus + H2O` has `κ = k [H2O] / N_A`. With the
macroscopic propensity, the mean encounter time of a pair at spacing `h` is
`(1 + h*/h)` times the macroscopic one. Below `h*` no propensity can match it.

**Mesh correction.** `[ssa] rate_correction = "hellander"` (the default) scales every
bimolecular propensity of a grid level of spacing `h > h*` by

```text
1 / (1 − h*/h)
```

This is eq. 6 of Hellander and Hellander (J. Chem. Phys. 152, 034104, 2020) with
`1/k = 1/k_r + 1/(4π σ D)`. It needs no reaction radius, and the mean encounter time
becomes the macroscopic one. Each coarsening level uses its own spacing. The factor is
applied when the propensity model of a level is built, so it costs nothing at run time.

- A reaction with `h ≤ h*` keeps the macroscopic propensity on that level. radkin warns
  once per such reaction on the fine grid. For diffusion-limited radical reactions `h*`
  is 1 to 4 nm. Reactions written with H2O have a large `k [H2O]`, so their `h*` is
  large: `e_aq + O_minus + H2O` gives 75 nm with the bundled chemistry. These stay
  uncorrected on nm grids.
- Reactions with a non-diffusing pair, first-order reactions and Michaelis–Menten
  reactions are not corrected.
- The formula assumes cubic voxels. An anisotropic grid runs uncorrected, with a
  warning, when the key is unset. An explicit `"hellander"` on an anisotropic grid is an
  error, and so is an explicit `"hellander"` on a unit-cell run, which has no grid.
- `rate_correction = "none"` gives the plain RDME. `snapshots.json` records which model
  ran.

Validation (`ssa::rdme::tests::mesh_correction_matches_the_association_time`): one
immobile A in the centre voxel and one B with a uniform start, `A + B → C`, on an
11³ grid with `h = 2 h*`. Over 4000 trajectories, the mean reaction time divided by
`V/κ` is 0.966 ± 0.015 with the correction and 1.458 ± 0.023 without; theory gives 1
and 1.5. The test runs 2000 trajectories. The few per cent left come from the finite,
zero-flux grid; Theorem 1 of Hellander et al. is for a periodic domain as N → ∞.

On `data/voxel_example.toml` (6 nm voxels, a synthetic 6-spur track, 100 ns; 64
trajectories, seed 1):

| G at 100 ns [/100 eV] | RDME, corrected (± SEM) | RDME, `"none"` (± SEM) | Deterministic |
|---|---|---|---|
| e_aq | 0.92 ± 0.05 | 1.13 ± 0.06 | 0.84 |
| OH_r | 2.18 ± 0.08 | 2.53 ± 0.07 | 2.10 |
| H3O_plus | 2.06 ± 0.06 | 2.57 ± 0.06 | 2.02 |
| H2O2 | 1.75 ± 0.04 | 1.77 ± 0.04 | 2.01 |
| H2 | 0.64 ± 0.04 | 0.68 ± 0.05 | 0.80 |

Most of the radical excess of the plain RDME came from this bias. With the correction,
e_aq, OH_r and H3O_plus are within 2 SEM of the deterministic run. H2O2 and H2, the
products of `2 A` recombinations, stay below it. With about one molecule per voxel, the
rate equations let a fraction of a molecule react with itself (`k [A]²` is not zero
for one molecule, `n (n − 1)` is). In this example four reactions with H2O stay
uncorrected (`h*` of 12 to 75 nm).

### Species levels

Hops dominate the cost of the RDME, and global coarsening keeps the whole grid fine as
long as any species is clustered. A numerous, fast background species, such as O2 in
aerated water (about 150 molecules per 100 nm voxel), then hops on the fine grid of
the radicals. `[ssa.species_levels]` puts a species on a coarser level of the same
2^L hierarchy (multi-grid RDME, ADR 0005; Cao and Erban 2014, Erban and Winkelmann
2024):

```toml
[ssa.species_levels]
O2 = 2      # O2 lives on 4x4x4 blocks of fine voxels; unset species stay on 0
```

- **Levels.** A species of level `L` holds its counts on the voxels of level `L` and
  hops there with `D / h_L²` per molecule and face. The hierarchy has as many levels as
  the deeper of `voxel.max_coarsen_levels` and the species levels; the grid dimensions
  must be divisible by `2^depth`. The key needs a voxel run with `solver = "ssa"`: unit
  cells and deterministic voxel runs refuse it. With every species on level 0 the run is
  the plain RDME, draw for draw.
- **Reactions.** A reaction belongs to the level of its coarsest reactant and fires in a
  voxel of that level, with that level's propensity, mesh correction included. The
  count of a finer reactant is summed over the fine voxels the coarse voxel contains:
  the pair reacts as if both occupied the coarse voxel (Hellander and Hellander 2020,
  III.C.4). A consumed molecule of the finer reactant is taken from one of the
  contained voxels, in proportion to their counts.
- **Products.** A product goes where the finest consumed molecule was: into the voxel
  of its own level that contains that spot or, if the product lives on a finer level
  than that molecule, uniformly into one of its voxels inside it (Erban and Winkelmann
  2024, eq. 20). For `e_aq + O2 → O2_minus` with O2 on level 2, O2_minus appears in the
  electron's voxel. Without a consumed molecule the reaction's voxel is the spot.
- **Global coarsening** is a floor: a species lives on the larger of its own level and
  the current global level. The dispersion test reads only the species that would merge
  (those on the current global level), and a seed refines the grid back to every
  species' own level. Initial solutions are drawn per voxel of each species' level;
  seeds land in the voxel of the species' level that contains the molecule; snapshots
  spread a coarse voxel evenly over its fine voxels. `snapshots.json` records the
  levels.
- **Tau-leaping** works with species levels (`src/ssa/rdme_multigrid_tau.rs`): the
  populations are the counts of each species on its level, cross-level reactions fire
  one firing at a time inside a leap with the rules above.
- **Suggestion.** For every background species (non-zero initial concentration, not
  seeded) the run logs `floor(log2 sqrt(D / D_slowest seeded))`, the level at which its
  hops per molecule match those of the slowest seeded species. It is a hint and is never
  applied.

Validation (`tests/ssa_rdme_multigrid.rs`, `tests/ssa_rdme_aerated.rs`, fixed seeds):

| Case | Compared | Result |
|---|---|---|
| A on level 1 converting into B on level 0 (8x2x2 grid) | mean of the 4 A and 32 B counts against the mean equations of the mgRDME (eq. 21-22), RK4 | worst 2.2 SEM, 400 runs |
| One e_aq on a 2 nm grid, one OH_r on level 2 (a single 8 nm voxel) | mean association time against `1 / c` of the level-2 voxel, mesh-corrected | 0.2 SEM, 4000 runs; the uncorrected `1 / c` is more than 10 SEM away |
| `large_count_case`, OH_r on level 1, `tau_epsilon = 0.03` | tau-leaping against the exact multi-grid run | 0.14 % of the species peak, 16 runs |
| Homogenising pure-water case, H2 on level 1, `max_coarsen_levels = 3` | species totals against the fine grid | worst 1.6 SEM, 24 runs |

Measured on an aerated track (pure water, 250 µmol/L O2, 8x8x8 grid of 100 nm voxels,
2 µs, 96 trajectories; `cargo test --release --test ssa_rdme_aerated -- --ignored
--nocapture`), against O2 on level 0, as the worst difference of the species totals over
about 1700 checks (15 species, 113 record times):

| Track | O2 level | Hops per trajectory | Fewer hops | Wall time | Worst difference |
|---|---|---|---|---|---|
| sparse: 300 each of e_aq, OH, H3O+, H in a 540 nm cube | 0 | 200 270 | | 0.93 s | |
| sparse | 1 | 48 584 | 4.1x | 0.30 s | 3.4 SEM |
| sparse | 2 | 14 121 | 14.2x | 0.14 s | 2.2 SEM |
| dense: 1500 each in a 140 nm cube | 0 | 205 682 | | 0.87 s | |
| dense | 1 | 54 642 | 3.8x | 0.42 s | 22 SEM |
| dense | 2 | 20 412 | 10.1x | 0.24 s | 28 SEM |

On the sparse track a few radicals per voxel meet 150 O2: O2 stays a background, and
the worst differences are those of noise over 1700 checks. On the dense track the
radicals outnumber the O2 of their voxels. The fine grid depletes O2 locally, a coarse
O2 voxel spreads the depletion over 8 or 64 times the volume, and e_aq and H react
with O2 too fast (HO2 6 molecules high at 20 ns). The level is therefore a choice of
accuracy against cost: a coarse species must stay close to uniform at the scale of its
voxels throughout the run. A track species is the worst case: OH- on level 1 in the
homogenising case leaves 40 % too much H3O+ at 1 µs, because its near
diffusion-limited recombination with H3O+ is spread over 8 times the volume.

## Cost

The cost is proportional to the number of events. On one core the exact SSA fires about
5 million events per second on the pure-water chemistry (73 reactions and 6 sources) and 15 million
on a three-reaction chemistry. The RDME fires about 3 million events per second, most of
them hops: a molecule hops `6 D / h²` times per second (8e8 per second for e_aq at
6 nm), so the cost grows as the number of molecules over `h²`. The fast acid-base equilibria fire continuously, so the
event count grows with the simulated time even when the slow chemistry is over.
