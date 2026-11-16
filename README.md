<p align="center">
  <img src="docs/assets/logo.svg" alt="radkin" width="520">
</p>

<p align="center">
  <a href="https://github.com/rtonneau/radkin/actions/workflows/ci.yml"><img src="https://github.com/rtonneau/radkin/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-LGPL--3.0--or--later-blue.svg" alt="License: LGPL-3.0-or-later"></a>
  <img src="https://img.shields.io/badge/rust-2024_edition-orange.svg?logo=rust" alt="Rust 2024 edition">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Windows-lightgrey.svg" alt="Platforms: Linux, Windows">
  <img src="https://img.shields.io/badge/C%20dependencies-none-success.svg" alt="No C dependency">
  <a href="https://github.com/rtonneau/radiopyo"><img src="https://img.shields.io/badge/validated%20against-radiopyo-6f42c1.svg" alt="Validated against radiopyo"></a>
</p>

<p align="center">
  <b>Chemical kinetics of water radiolysis, from a single track to the homogeneous solution.</b>
</p>

---

radkin simulates the chemistry of irradiated water under conventional and
ultra-high-dose-rate (FLASH) irradiation. It is a Rust port of
[radiopyo](https://github.com/rtonneau/radiopyo) and reproduces its results, then extends
it to 3D voxel geometries seeded by Geant4 tracks and to stochastic chemistry.

- **Chemistry in TOML.** Acid/base, equilibrium, Arrhenius, constant and
  Michaelis–Menten reactions; radiopyo configuration files work unchanged.
- **Radiolytic yields.** G-values from the Boscolo 2020 LET fit, with an optional
  dose-rate factor.
- **Beams.** Constant, single-pulse, pulsed and macro-pulsed.
- **Geant4 coupling.** Geant4 provides the track up to the end of the physico-chemical
  stage; radkin runs the whole chemical stage from there, in place of Geant4-DNA's
  chemistry.
- **Reproducible and parallel.** Seeded ChaCha8 streams give the same output on every
  platform and thread count; voxels and stochastic trajectories run in parallel (rayon).

## Three solvers

| | Solver | Model | Use it for |
|---|---|---|---|
| 1 | **LSODA** (deterministic) | Rate equations on a well-mixed cell or a voxel grid (Strang splitting with ADI diffusion), solved by [rsoda](https://github.com/rtonneau/rsoda) | Homogeneous kinetics, dose-rate studies, long times |
| 2 | **Gillespie SSA / RDME** (stochastic) | Integer molecule counts; exact SSA or tau-leaping in a cell, reaction-diffusion master equation on a grid | Small volumes, low doses, a few molecules per voxel |
| 3 | **IRT** (particle-based) | Every molecule is a particle; independent reaction times per pair, with Geant4-DNA's formulas | Track structure from 1 ps, before diffusion has spread the spurs |

The solvers chain: an IRT stage can run the first microsecond of a track and hand its
molecules over to the RDME grid.

### Gillespie algorithm and RDME features

- **Exact SSA** (Gillespie's direct method) and **tau-leaping** (Cao, Gillespie and
  Petzold 2006) on the unit cell, with propensities built from the same matrices as the
  LSODA kernel.
- **Reaction-diffusion master equation** on the voxel grid with the next-subvolume
  method: reactions fire inside voxels, diffusion is a hop to a neighbour at rate `D/h²`.
- **Mesh correction** of bimolecular rates (Hellander, Hellander and Petzold 2012), so
  that nanometre voxels keep the macroscopic rate; reversible pairs share one factor and
  keep their equilibrium constant.
- **Adaptive coarsening**: the grid merges 2×2×2 blocks once the counts are
  statistically uniform and refines back before each new track.
- **Multi-grid species levels**: numerous background species such as O₂ live on a
  coarser grid, which cuts the number of hops by an order of magnitude.
- **Tau-leaping on the grid** for reactions and hops together.
- **Background run**: a seedless ensemble is subtracted so that G-values are net of the
  background chemistry.
- **Ensembles** of independent trajectories, reported as mean and standard deviation.

## Quick start

```bash
cargo build --release          # target/release/radkin
cargo test --release           # unit tests and golden comparisons with radiopyo
```

```bash
# Homogeneous kinetics (LSODA)
radkin run data/config_Labarbe.toml --t-span 1e-7 10 --out out/labarbe

# Stochastic (Gillespie SSA); the run file needs an [ssa] section with its volume
radkin run my_ssa.toml --t-span 0 1e-3 --solver ssa --out out/ssa

# 3D voxel grid seeded by a Geant4 track (LSODA, or --solver ssa for the RDME)
radkin voxel data/voxel_example.toml --out out/voxel

# Particle-based IRT stage
radkin irt data/irt_example.toml --out out/irt
```

As a library:

```rust
use radkin::unit_cell::{RunOptions, UnitCell};

let mut cell = UnitCell::from_toml("data/config_Labarbe.toml")?;
let res = cell.run((1e-7, 10.0), &RunOptions::default(), None)?;
let (labels, rows) = res.concentrations();
```

Results are written as CSV, JSON and NumPy (`.npy`, `.npz`) files. Logging:
`RUST_LOG=info`.

## Validation

| Against | Result |
|---|---|
| radiopyo | Right-hand side, Jacobian and G-values at 1e-12; 27-case acceptance matrix (3 chemistries × 4 beams × run modes) within `rtol = 1e-6`, observed differences ~1e-8 |
| LSODA | SSA ensemble means and RDME totals agree with the deterministic runs at large counts |
| Geant4-DNA IRT_syn | On 2000 single-track events, G(t) within 1 % for e⁻aq, •OH, H•, H₃O⁺, OH⁻ and H₂ from 100 ps to 10 µs, within 4 % for H₂O₂ |

## Documentation

| Topic | |
|---|---|
| Deterministic solver: command line, outputs, library, comparison with radiopyo | [docs/lsoda.md](docs/lsoda.md) |
| Voxel reaction-diffusion: configuration, coarsening, checkpoints, caveats | [docs/voxel.md](docs/voxel.md) |
| Stochastic solver and RDME: usage, model, validation, cost | [docs/ssa.md](docs/ssa.md) |
| IRT stage: model, pair rules, boundaries, handover, validation | [docs/irt.md](docs/irt.md) |
| Performance measurements | [docs/perf/](docs/perf/) |
| Bundled data | [data/README.md](data/README.md) |

## Reference chemistry

`data/tonneau_2025.toml` is the chemistry of Tonneau et al. 2025, *Phys. Med. Biol.*
70 235021, with the inputs of the paper's runs. `data/chem_g4.toml` is the Geant4-DNA
option3 chemistry used by the IRT stage.

## License

[LGPL-3.0-or-later](LICENSE), as radiopyo.
