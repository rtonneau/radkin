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

## Cost

The cost is proportional to the number of events. On one core the exact SSA fires about
5 million events per second on the pure-water chemistry (73 reactions and 6 sources) and 15 million
on a three-reaction chemistry. The fast acid-base equilibria fire continuously, so the
event count grows with the simulated time even when the slow chemistry is over.
