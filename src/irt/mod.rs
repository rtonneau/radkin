//! IRT stage: every molecule is a particle with a position, and reactions are
//! drawn as independent reaction times per pair, and per molecule for
//! background reactions (ADR 0007, `docs/irt.md`).

pub mod model;
pub mod sample;
pub mod species;

pub use model::{IrtModel, MoleculeReaction, PairReaction, TypeSource};
pub use species::{SpeciesProps, SpeciesTable};
