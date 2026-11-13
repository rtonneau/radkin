//! radkin: Rust port of radiopyo (homogeneous water radiolysis kinetics),
//! integrated with the LSODA solver of `rsoda`.

pub mod chemistry;
pub mod config;
pub mod consts;
pub mod error;
pub mod physics;
pub mod poly;
pub mod pyfloat;

pub use error::{RadkinError, Result};
