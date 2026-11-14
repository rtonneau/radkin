//! Error type of the crate.

use std::path::PathBuf;

/// Every error radkin reports.
#[derive(Debug, thiserror::Error)]
pub enum RadkinError {
    /// A reaction or species string could not be parsed.
    #[error("unable to parse '{input}': {msg}")]
    Parse { input: String, msg: String },
    /// Invalid or inconsistent configuration.
    #[error("CONFIG::{0}")]
    Config(String),
    /// File access error.
    #[error("{path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// Invalid beam definition.
    #[error("BEAM::{0}")]
    Beam(String),
    /// A configured beam ON time is shorter than `MIN_BEAM_ON_TIME`.
    #[error("ShortBeamOnError: {0}")]
    ShortBeamOn(String),
    /// A run was requested without any beam.
    #[error("No beam is defined for the simulation")]
    NoBeamDefined,
    /// The LSODA solver failed on a span.
    #[error(
        "LSODA failed on t_span=({start}, {stop}): last reached t={t_reached:.6e} s; rsoda message: {message}"
    )]
    SolverFailed {
        start: f64,
        stop: f64,
        t_reached: f64,
        message: String,
    },
    /// G-value lookup failure (e.g. LET out of range).
    #[error("G-values: {0}")]
    GValue(String),
    /// Unreadable or incompatible voxel checkpoint.
    #[error("CHECKPOINT::{0}")]
    Checkpoint(String),
    /// Invalid value given to a function.
    #[error("{0}")]
    Value(String),
}

/// Result alias of the crate.
pub type Result<T> = std::result::Result<T, RadkinError>;

impl RadkinError {
    /// Io error for `path`.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }
}
