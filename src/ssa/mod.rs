//! Stochastic solver (`solver = "ssa"`): the reaction network on integer
//! molecule counts, simulated with Gillespie's algorithm. The exact SSA is
//! the reference; tau-leaping is an opt-in approximation of it. Both take
//! the same chemistry as the LSODA path.

pub mod config;
pub mod ensemble;
pub mod exact;
pub mod model;
pub mod output;
pub mod rdme;
pub mod rdme_coarsen;
pub mod rdme_multigrid;
pub mod rdme_multigrid_tau;
pub mod rdme_tau;
pub mod rng;
pub mod schedule;
pub mod tau;
pub mod voxel_output;
pub mod voxel_sim;

pub use config::{RateCorrection, SolverConfig, SolverKind, SsaConfig};
pub use ensemble::{EnsembleOptions, SsaResult, run_ensemble};
pub use exact::{Cell, Exact, Segment, Stepper, Trajectory, initial_counts, run_trajectory};
pub use model::{Law, Reaction, Slot, Source, SsaModel};
pub use output::{SsaRunMeta, write_ssa_dir};
pub use rdme::{Rdme, RdmeSetup, ReversiblePair};
pub use rdme_coarsen::dispersion_threshold;
pub use rdme_multigrid::{Multigrid, MultigridSetup};
pub use rng::{SsaRng, det_ln, ln_factorial};
pub use schedule::{beam_segments, delivered_dose};
pub use tau::TauLeap;
pub use voxel_output::{rdme_yields_csv, write_rdme_dir};
pub use voxel_sim::{BackgroundRun, RdmeOutput, run_rdme};
