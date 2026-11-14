//! Stochastic solver (`solver = "ssa"`): the reaction network on integer
//! molecule counts, simulated with Gillespie's algorithm. The exact SSA is
//! the reference; it takes the same chemistry as the LSODA path.

pub mod config;
pub mod ensemble;
pub mod exact;
pub mod model;
pub mod output;
pub mod rng;
pub mod schedule;

pub use config::{SolverConfig, SolverKind, SsaConfig};
pub use ensemble::{EnsembleOptions, SsaResult, run_ensemble};
pub use exact::{Cell, Exact, Segment, Stepper, Trajectory, initial_counts, run_trajectory};
pub use model::{Law, Reaction, Slot, Source, SsaModel};
pub use output::{SsaRunMeta, write_ssa_dir};
pub use rng::{SsaRng, det_ln};
pub use schedule::{beam_segments, delivered_dose};
