//! Constants shared with radiopyo (same values).

use std::path::PathBuf;

/// Beginning of every beam ON phase skipped, its radiolytic yield added analytically [s].
pub const T_ZERO: f64 = 1e-7;
/// Default first step given to LSODA [s].
pub const FIRST_STEP: f64 = 1e-15;
/// Shortest configured beam ON time [s] (0.4 µs; a shorter configured ON time is rejected, a shorter run-time ON span is added analytically).
pub const MIN_BEAM_ON_TIME: f64 = 0.4e-6;
/// Relative slack on the configured ON time check.
pub const ON_TIME_RTOL: f64 = 1e-9;
/// Relative tolerance used when counting periods.
pub const COUNT_RTOL: f64 = 1e-9;
/// Half-saturation concentration of the depletion switch [mol/L].
pub const DEPLETION_LIMIT: f64 = 1e-12;
/// Default relative tolerance of the solver.
pub const DEFAULT_RTOL: f64 = 1e-10;
/// Default absolute tolerance of the solver [µmol/L].
pub const DEFAULT_ATOL: f64 = 1e-12;
/// Default transition time after the max dose in single-pulse runs [s].
pub const DEFAULT_TRANS_TIME: f64 = 1e-4;
/// First thinned point offset, as a fraction of the chunk span.
pub const THINNING_FIRST_OFFSET: f64 = 1e-6;
/// Molar gas constant [J/mol/K] (scipy.constants.R).
pub const GAS_CONSTANT: f64 = 8.314462618;
/// Elementary charge [C].
pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;
/// Avogadro constant [1/mol].
pub const AVOGADRO: f64 = 6.02214076e23;
/// Default temperature [degC].
pub const DEFAULT_TEMPERATURE: f64 = 25.0;
/// 0 degC in kelvin.
pub const ZERO_CELSIUS: f64 = 273.15;
/// Concentration of water [mol/L] as a polynomial of the temperature [degC].
pub const H2O_CC: [f64; 5] = [55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9];
/// Fit of the dose-rate dependence of the G-values: `GF_A * dr^GF_B`.
pub const GF_A: f64 = 0.97207704;
/// Exponent of the dose-rate G-value fit.
pub const GF_B: f64 = -0.10972;

/// Folder of the bundled data files (`RADKIN_DATA_DIR` overrides it).
pub fn data_dir() -> PathBuf {
    match std::env::var_os("RADKIN_DATA_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data"),
    }
}
