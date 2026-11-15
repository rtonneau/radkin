//! IRT stage: every molecule is a particle with a position, and reactions are
//! drawn as independent reaction times per pair, and per molecule for
//! background reactions (ADR 0007, `docs/irt.md`).

pub mod config;
pub mod model;
pub mod output;
pub mod sample;
pub mod search;
pub mod sim;
pub mod solver;
pub mod species;

pub use config::{BoundaryKind, IrtConfig, Resampling};
pub use model::{IrtModel, MoleculeReaction, PairReaction, TypeSource};
pub use output::{irt_yields_csv, write_irt_dir};
pub use sim::{IrtOutput, IrtSeed, IrtSimulation, run_irt, run_irt_with};
pub use solver::{Boundary, Stage, StageOptions, StageStats};
pub use species::{SpeciesProps, SpeciesTable};
