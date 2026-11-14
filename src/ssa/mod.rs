//! Stochastic solver (`solver = "ssa"`): the reaction network on integer
//! molecule counts, simulated with Gillespie's algorithm. The exact SSA is
//! the reference; it takes the same chemistry as the LSODA path.

pub mod exact;
pub mod model;
pub mod rng;
pub mod schedule;

pub use exact::{Cell, Exact, Segment, Stepper, Trajectory, initial_counts, run_trajectory};
pub use model::{Law, Reaction, Slot, Source, SsaModel};
pub use rng::{SsaRng, det_ln};
pub use schedule::{beam_segments, delivered_dose};
