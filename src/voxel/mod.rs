//! Phase 2: voxelised geometry with species diffusion, coupled to the
//! per-cell chemistry by Strang splitting (ADR 0002).

pub mod config;
pub mod diffusion;
pub mod grid;
pub mod reaction;

pub use config::{DEFAULT_DIFFUSION, VoxelConfig};
pub use diffusion::{Diffusion, Scheme};
pub use grid::{VoxelGrid, VoxelState};
