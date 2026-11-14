//! Phase 2: voxelised geometry with species diffusion, coupled to the
//! per-cell chemistry by Strang splitting (ADR 0002).

pub mod checkpoint;
pub mod config;
pub mod diffusion;
pub mod driver;
pub mod grid;
pub mod output;
pub mod reaction;
pub mod seed;
pub mod sim;

pub use checkpoint::{CHECKPOINT_VERSION, Checkpoint, SeedRecord};
pub use config::{DEFAULT_DIFFUSION, VoxelConfig};
pub use diffusion::{Diffusion, Scheme};
pub use driver::{StepControl, StepRecord, VoxelRun};
pub use grid::{VoxelGrid, VoxelState};
pub use output::write_voxel_dir;
pub use seed::{
    DEFAULT_GEANT4_ALIASES, LoadedSeed, SeedPolicy, SeedReport, SeedSpec, load_seeds, run_schedule,
};
pub use sim::{RunControl, RunSettings, Snapshot, VoxelOutput, VoxelSimulation};
