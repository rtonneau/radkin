//! radkin: Rust port of radiopyo (homogeneous water radiolysis kinetics),
//! integrated with the LSODA solver of `rsoda`.

pub mod chemistry;
pub mod config;
pub mod consts;
pub mod error;
pub mod gvalues;
pub mod kernel;
pub mod matrices;
pub mod physics;
pub mod poly;
pub mod pyfloat;
pub mod solver;

pub use error::{RadkinError, Result};
