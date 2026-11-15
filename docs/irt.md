# IRT stage

In the IRT stage every molecule is a particle with a position, and reactions are drawn
as independent reaction times: one per pair of particles, and one per particle for the
reactions it has alone or with the background. It runs from the molecules of a track at
1 ps to a handover time, where the molecules are binned onto the voxel grid, or to the
end of the run. The code is in `src/irt/`; the decisions are in ADR 0007
(`.work/adr/0007-irt-stage.md`). This file records how the stage is built and how it
was validated.

```bash
radkin irt data/irt_example.toml --out out/irt        # the IRT stage alone
radkin voxel CONFIG --solver ssa --out out/run        # IRT stage, then the RDME (with [irt])
```

## Why

- A study with Geant4-DNA needed one Geant4 pass per handover time, and Geant4's IRT
  runs one track at a time, without background reactions.
- Tracks overlap before any practical handover above about 1e7 Gy/s.
- The RDME does not reproduce the chemistry of a spur from 1 ps, whatever the voxel
  size: a pair 1 to 3 nm apart is well mixed inside a voxel (PR #18: •OH +40 % at
  10 us on single tracks).

## Model of a network

`IrtModel` (`src/irt/model.rs`) sorts the mass-action reactions of a network by their
dynamic reactants:

| Reactants | In the IRT stage |
|---|---|
| two dynamic species | a **pair reaction** between two particles |
| one dynamic species (first order, or with fixed species such as water) | a **molecule reaction** of rate `k · Π[fixed]` |
| a pair reaction whose partner has an initial concentration | also a molecule reaction with that **background species**, of rate `k c` (`2 k_obs c` for `A + A`) |
| no dynamic species (water autoionisation) | skipped: the background keeps itself at equilibrium |

More than two dynamic reactants, "negative" products and Michaelis–Menten reactions are
refused. A configuration with a `[beam]` is refused too: the molecules of the stage come
from seeds, not from G-values.

`k_obs` is the rate constant of the reaction events (`events = k_obs [A][B]`, and
`k_obs [A]²` for a self-reaction). It is the value radkin stores, and Geant4's observed
rate constant. A scalar `k_value` of `2 A -> …` in a TOML is divided by 2 when it is
read, so the Geant4 export writes `2 k_obs` there.

### Reaction types and radii

The formulas are those of Geant4-DNA 11.4.3 (`G4DNAMolecularReactionTable.cc`).

| | Fully diffusion-controlled | Partially diffusion-controlled |
|---|---|---|
| reaction radius `R` | the effective radius | the sum of the two van der Waals radii |
| effective radius `r_eff` | `k_obs / (4π D N_A)`, `D` the sum of the two coefficients (the single one for `A + A`) | `R`, or `−r_c / (1 − exp(r_c / R))` for two ions, `r_c` the Onsager radius (293.15 K, permittivity 80.1) |
| `k_dif` | `k_obs` | `4π D r_eff N_A`, halved for `A + A` |
| `k_act` | infinite | `k_dif k_obs / (k_dif − k_obs)` |
| contact probability | 1 | `R_s / (R_s + (k_dif / k_act)(r_eff + R_s))`, `R_s` = 0.29 nm |

The type of a reaction is:

1. the `irt_type` key of its entry (`"full"` or `"partial"`), which
   `tools/export_g4_chem.py` writes from Geant4's reaction type;
2. else `"full"` if `k_obs >= k_dif` and `"partial"` otherwise, with `k_dif` from the
   species table;
3. `"full"`, with a warning, if a reactant is not in the species table.

On the 52 pair reactions of Geant4's option3 the rule of point 2 gives Geant4's type for
47. The other 5 (`H + H`, `H + e_aq`, `H + O`, `H + O⁻`, `e_aq + e_aq`) are typed fully
diffusion-controlled by Geant4 although their rate constant is below the diffusion
limit of their radii. Geant4 never types "partial" a reaction at or above its limit.
With the `irt_type` keys, the 52 effective and reaction radii equal the values Geant4
dumped at run time to 1e-9 (`tests/irt_model.rs`).

### Species table

`SpeciesTable::water()` (`src/irt/species.rs`) holds the radius, the charge and the
diffusion coefficient of the 15 water-radiolysis species, with Geant4's values. A
configuration can override or add entries (`[irt.species.<label>]`); the
`[voxel.diffusion]` coefficients of a voxel run apply to the stage too.

## Sampling kernel

`src/irt/sample.rs`, from `G4DiffusionControlledReactionModel::GetTimeToEncounter` and
`G4DNAMakeReaction`:

- **Reaction time of a pair** at the distance `r`. The pair reaches `r_eff` with the
  probability `W = r_eff / r`, at the time `(r − r_eff)² / (4 D erfcinv(U / W)²)`. A
  partially diffusion-controlled pair then reacts with the probability
  `k_act / (k_act + k_dif)`, after a second time. Together they give exactly the law of
  a radiation boundary, not only its long-time limit.
- **Contact rule.** A pair injected inside its effective radius reacts at once with the
  contact probability of its reaction (`start_contact`).
- **Positions.** The reactants are put `r_eff` apart; the angle of their separation
  with the initial one, and the diffusion of the pair over the step, are sampled. Every
  product of a reaction gets the same position, between the reactants.

Two points differ from Geant4's code, where it departs from the formulas it implements
(Plante and Devroye 2017): the pair diffuses with the standard deviation `sqrt(σ)`
(Geant4 passes the variance, which in its units is about zero), and the angle is taken
from the initial separation (Geant4 takes it from the z axis).

The kernel does not allocate and uses IEEE basic operations, the square root and
radkin's own logarithm, exponential and normal quantile (`src/ssa/rng.rs`), so a seed
gives the same draws on every platform. Checked against the analytic reaction
probabilities at 5 to 6 times from 1 ps to 1 us, for fully and partially
diffusion-controlled pairs and for ion pairs (`tests/irt_sample.rs`, 400 000 pairs,
4 standard errors).

## The two pair rules

### Full resampling (default)

`src/irt/solver.rs`. One pass draws a reaction time for every pair inside the cutoff
radius, and one time for the reactions of single particles, and keeps the earliest. The
event moves the time, puts the reactants and products at sampled positions, moves every
other particle by its free diffusion over the step, and discards every time: the next
pass draws them all again.

- The cutoff radius is `r_max + 2 sqrt(D_max t_left) erfcinv(p)`, with `p` =
  `cutoff_probability` (0.01, as Geant4) and `t_left` the time to the end of the stage.
  Only particles of a species that has a partner among the particles are searched.
- Record times are not stops: particles move at reactions, and a record reads the
  counts.
- A pair that a move brought inside its effective radius gets no reaction time in that
  pass: its sampled time said it had not reacted. It bounds the next step to
  `r_eff² / D`, after which every particle is moved and looked at again.

### Incremental resampling

`src/irt/incremental.rs`, `resampling = "incremental"`. Every reaction time is kept in
a queue. A time is drawn for a pair once, when the later of its two particles appears
(injected or created), and for a particle alone when it appears. After a reaction only
the pairs of its products are drawn.

- A particle keeps the position of the last time it was looked at. It is moved to the
  present when a new particle looks for partners around it, and when it reacts alone.
- A pair that fires is put at its reaction positions from the two positions its time
  was drawn from.
- Between reflective walls every pair is drawn again at the look-ahead limit of the box
  (below), as the full rule does.

The two rules agree on the same seeds, for every species at every record within 4
standard errors: single tracks, five overlapping tracks with the pH 7 background, and a
box with reflective walls (`tests/irt_incremental.rs`). Geant4's IRT_syn keeps its
reaction set between reaction steps and clears it only after a step limited by a
process, so it is closer to the incremental rule.

## Background reactions

A species with an initial concentration (O2 in aerated water, H3O⁺ and OH⁻ at the pH of
the configuration) is a background species. A particle reacts with it at the rate
`k c`, with an exponential reaction time, like a first-order reaction or a reaction
with water.

- The background keeps its concentration. The molecules it loses are counted
  (`background_consumed` in `irt.json`), not debited.
- A background species is not made of particles; the molecules of the same species a
  track makes are particles.
- `background_reactions = false` turns off every reaction of one particle.

Checked: a first-order decay and the scavenging of e_aq by 250 uM of O2 follow
`exp(−k t)`; a self-reaction with its own background has the rate `2 k_obs c`
(`tests/irt_background.rs`).

## Boundaries

`boundary` and a box (`[irt.box]`, or the grid of `[voxel.grid]`):

- **`"none"`**: free space.
- **`"reflective"`**: every position is folded back into the box, and a pair near a
  wall also samples a reaction time against the mirror image of its partner, across
  every wall within the cutoff radius, with the edge and corner images. Only the first
  image across each wall is used, so a pass looks ahead no further than the time for
  which its cutoff radius is the shortest edge of the box (22 us for a 2 um box); past
  it every particle is moved and every pair is drawn again.
- **`"outflow"`**: a molecule that a move puts outside the box is removed. This is what
  Geant4 does at the faces of its water volume, and it exists to compare with a Geant4
  run made in a small box. It is not an absorbing wall: about half as many molecules
  are lost when the moves are long.

Checked (`tests/irt_walls.rs`):

- a uniform distribution stays uniform between the walls, and molecules started in a
  corner fill the box;
- one fully diffusion-controlled pair 4 nm (7 reaction radii) from a wall of a 60 nm
  box: survival at 20 ns 0.882 with the image, 0.889 by a Brownian walk with
  reflection, 0.914 without the image (5 standard errors off);
- 30 + 30 molecules of a partially diffusion-controlled `A + B` in a 40 nm box follow
  the exact stochastic rate equation.

**Limit of the image rule.** In a box of 11 reaction radii the walk lies between the
image rule (reaction rate 1.5 times too high) and no image (0.85 times). A sink and its
mirror image compete for the same molecules, which independent pairs cannot see: near a
wall at the distance `h` the true rate is about `1 / (1 + r_eff / 2h)` times the free
one, while the image rule restores the free rate. The rule is for walls far compared
with the reaction radius (0.5 nm), which is every box radkin is meant for.

## Seeds and scheduled tracks

`[[irt.seed]]` takes the keys of `[[voxel.seed]]` (`file`, `time`, `events`,
`offset_nm`, `energy_ev`). A track is injected at its time, after every particle has
been moved to that time, and its molecules join the pair search. `G` in `yields.csv`
divides by the energy injected so far. Tracks 50 um apart and 0.2 us one after the
other reproduce the sum of the single-track counts at their ages
(`tests/irt_schedule.rs`).

`start_contact` (default `true`) applies the contact rule to injected molecules. Set it
to `false` for molecules Geant4 wrote after its chemistry started (its seed pass stops
at 1.001 ps, after the contact rule at 1 ps): on such seeds `true` gives too much H2O2
at 10 ps.

## Handover to the RDME

`radkin voxel` with an `[irt]` section and the ssa solver runs, in each trajectory, the
IRT stage to `irt.t_end`, bins its molecules onto the fine grid and carries on with the
RDME (`run_rdme_from_irt`, `src/ssa/voxel_sim.rs`).

- The stage needs `boundary = "reflective"`, with the grid as its box.
- The seeds are those of `[[irt.seed]]`, or of `[[voxel.seed]]` when `[irt]` has none.
  A seed at or after the handover time is an error.
- The ensemble is that of `[ssa]` (`trajectories`, `seed`).
- The grid starts at the handover time, from the Poisson draw of the initial solution.
  `voxel.run.t_start` is not used; a snapshot before the handover is an error.
- `yields.csv` starts at the start of the stage. Before the handover `N` is the
  particles plus the grid's draw of the initial solution; `G` is net of it.
  `snapshots.json` gets an `irt_stage` object.
- No molecule is lost: the net molecules before the handover equal the counts of the
  same stage run on its own, and at the handover the grid holds every one
  (`tests/irt_handover.rs`, both pair rules).
- The deterministic voxel solver has no IRT stage.

## Configuration

```toml
[irt]
t_end = 1e-6             # s (required): end of the stage; the handover time in `radkin voxel`
t_start = 1e-12          # s (default: the earliest seed time)
trajectories = 16        # default 1 (`radkin irt` only)
random_seed = 1          # default 0 (`radkin irt` only)
resampling = "full"      # "full" or "incremental"
start_contact = true     # contact rule for injected molecules
background_reactions = true
boundary = "none"        # "none", "reflective" or "outflow" (the last two need a box)
wall_images = true       # mirror images at reflective walls
cutoff_probability = 0.01
record_times = [1e-9]    # s (default: 100 log-spaced times)
clip = false             # drop seed molecules outside the box instead of failing

[irt.box]                # default: the box of [voxel.grid]
lo_nm = [-1000, -1000, -1000]
hi_nm = [1000, 1000, 1000]

[[irt.seed]]             # default: the [[voxel.seed]] entries
file = "molecules.csv"
time = 1e-12
events = [0]
offset_nm = [0, 0, 0]
energy_ev = 615.0

[irt.geant4]             # as [voxel.geant4]
ignore = ["H2O^0"]

[irt.species.Cl_minus]   # overrides or adds an entry of the species table
radius_nm = 0.18
charge = -1
diffusion = 2.0e-9       # m^2/s
```

## Outputs of `radkin irt`

| file | content |
|---|---|
| `yields.csv` | `t`, `N_<sp>` (mean particles), `G_<sp>` when seeds give `energy_ev`, then `N_<sp>_std` and `G_<sp>_std` (standard deviation over the trajectories): the columns of the RDME `yields.csv` |
| `irt.json` | settings, event counts, firings of each reaction, background consumed, molecules that left the box, the type and effective radius of every pair reaction, warnings |

## Validation against Geant4-DNA IRT_syn

`python -I tools/compare_irt_geant4.py --events 100 [--boundary outflow]` runs each of
the 100 events of the reference (100 keV electrons through a 2 um water box, option4
physics, option3 chemistry, 615 eV per event; 9491 molecules at 1.001 ps) with 16
trajectories and compares the molecule counts at Geant4's 71 checkpoints. `z` is the
difference over its standard error (the chemistry noise of both codes and the
difference of their two initial states). Report of 2026-11-15:
`.scratch/tests/2026-11-15__irt-syn-port/03-report.md`.

- **Free space, to 300 ns:** RMS of z = 0.93 over 825 values; every value of the decade
  table within 3 standard errors.
- **The Geant4 reference is not free space after 300 ns.** It loses molecules at the
  faces of its box: its net charge drifts from −5 to −83 and its H2 decreases.
- **With the same box (`boundary = "outflow"`), to 10 us:** RMS of z = 1.00 over 1065
  values with the full rule and 1.05 with the incremental rule. Every value of the
  decade table is within 3 standard errors, except H3O⁺ at the last record with the
  full rule (−3.2).

G in molecules / 100 eV, radkin (full rule) / Geant4, with the outflow box:

| species | 10 ps | 1 ns | 100 ns | 1 us | 10 us |
|---|---|---|---|---|---|
| e_aq | 4.284 / 4.293 | 3.747 / 3.748 | 2.777 / 2.757 | 2.547 / 2.506 | 2.155 / 2.137 |
| •OH | 5.256 / 5.281 | 3.925 / 3.987 | 2.869 / 2.939 | 2.693 / 2.752 | 2.400 / 2.497 |
| H• | 1.081 / 1.061 | 0.892 / 0.917 | 0.897 / 0.930 | 0.884 / 0.925 | 0.750 / 0.793 |
| H2O2 | 0.091 / 0.075 | 0.485 / 0.463 | 0.661 / 0.629 | 0.667 / 0.637 | 0.634 / 0.627 |
| H2 | 0.127 / 0.130 | 0.193 / 0.197 | 0.324 / 0.323 | 0.354 / 0.348 | 0.338 / 0.338 |

In free space the same tracks give, at 10 us: e_aq 2.51, •OH 2.66, H• 0.91, H2O2 0.68,
H2 0.38. A Geant4 reference whose molecules cannot leave the water would allow the
comparison in free space over the last decade.

`tests/irt_single_track.rs` runs ten of these events (`tests/data/irt/`) in the test
suite.

## Cost

`docs/perf/irt-scaling.md`. To a 1 us handover, one trajectory, one thread:

| overlapping tracks | molecules | incremental | full |
|---|---|---|---|
| 10 | 869 | 6 ms | 0.33 s |
| 30 | 2 607 | 47 ms | 7.7 s |
| 100 | 8 690 | 0.47 s | 281 s |
| 1000 | 86 900 | 62 s | not run |

The 100 single-track events of the comparison take 6 s with the full rule and 1.6 s
with the incremental rule (16 trajectories each, 16 cores).

## Limits

- **Independent pairs.** IRT draws the fate of every pair on its own. The error grows
  as tracks overlap densely, and nothing here measures it. Results with several
  overlapping tracks, and aerated results, have no reference at the particle level yet:
  treat them as provisional until a Geant4 multi-track reference exists (idea
  `geant4-uhdr-multi-track-reference`).
- **Background.** The background is not consumed. A stage long or dense enough to
  deplete a scavenger needs an earlier handover to the grid.
- **Full resampling** draws a pair again as if no depletion had built up around it. When
  it is repeated many times before a reaction (a small box, with a move every
  look-ahead limit) it raises a diffusion-limited rate by about
  `r_eff / sqrt(π D Δt)`.
- **Networks that do not come from Geant4** run, with the type rule above, but have no
  reference at the particle level.
- **One handover time per run.** A pulse longer than the handover time needs a later
  handover, or a run of the IRT stage alone.
- The stage starts at 1 ps: the dissociation of water and the thermalisation of the
  electron stay in Geant4.
