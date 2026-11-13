//! Phase 2: voxelised geometry with species diffusion, coupled to the
//! per-cell chemistry by Strang splitting (ADR 0002).

pub mod config;
pub mod grid;

pub use config::{DEFAULT_DIFFUSION, VoxelConfig};
pub use grid::{VoxelGrid, VoxelState};
