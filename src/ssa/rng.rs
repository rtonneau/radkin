//! Random numbers of the stochastic solver.
//!
//! One ChaCha8 generator per trajectory: the config seed keys the generator
//! and the trajectory index selects its stream, so trajectories are
//! independent by construction and do not depend on the thread that runs
//! them. The samplers use IEEE basic operations only (the logarithm is
//! [`det_ln`], not the platform's libm), so a seed gives the same draws on
//! every platform.

use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};

/// 2^-53.
const ULP53: f64 = 1.0 / 9_007_199_254_740_992.0;

/// Natural logarithm of a finite `x > 0`, from IEEE basic operations only
/// (absolute error below 1e-15 * max(1, |ln x|)).
///
/// `x = m * 2^e` with `m` in `(sqrt(1/2), sqrt(2)]`, then
/// `ln m = 2 atanh(s)`, `s = (m - 1) / (m + 1)`, `|s| < 0.172`.
pub fn det_ln(x: f64) -> f64 {
    debug_assert!(x > 0.0 && x.is_finite(), "det_ln({x})");
    // Subnormals: scale into the normal range first.
    let (x, shift) = if x < f64::MIN_POSITIVE {
        (x * 18_014_398_509_481_984.0, -54) // 2^54
    } else {
        (x, 0)
    };
    let bits = x.to_bits();
    let mut e = ((bits >> 52) & 0x7ff) as i64 - 1023 + shift;
    let mut m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    if m > std::f64::consts::SQRT_2 {
        m *= 0.5;
        e += 1;
    }
    let s = (m - 1.0) / (m + 1.0);
    let s2 = s * s;
    // 1 + s^2/3 + s^4/5 + ... + s^26/27 (Horner); the next term is below 1e-22.
    let mut p = 1.0 / 27.0;
    let mut k = 25.0;
    while k >= 1.0 {
        p = p * s2 + 1.0 / k;
        k -= 2.0;
    }
    2.0 * s * p + e as f64 * std::f64::consts::LN_2
}

/// Generator of one trajectory.
#[derive(Debug, Clone)]
pub struct SsaRng(ChaCha8Rng);

impl SsaRng {
    /// Generator of trajectory `stream` for the config seed `seed`.
    pub fn new(seed: u64, stream: u64) -> Self {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        rng.set_stream(stream);
        Self(rng)
    }

    /// Uniform in `[0, 1)` (53 bits).
    pub fn uniform(&mut self) -> f64 {
        (self.0.next_u64() >> 11) as f64 * ULP53
    }

    /// Uniform in `(0, 1]`.
    pub fn uniform_open(&mut self) -> f64 {
        ((self.0.next_u64() >> 11) + 1) as f64 * ULP53
    }

    /// Exponential of mean 1.
    pub fn exp1(&mut self) -> f64 {
        -det_ln(self.uniform_open())
    }

    /// Uniform integer in `0..n` (`n >= 1`), by rejection (no modulo bias).
    pub fn below(&mut self, n: u64) -> u64 {
        debug_assert!(n >= 1);
        let zone = u64::MAX - (u64::MAX % n + 1) % n;
        loop {
            let v = self.0.next_u64();
            if v <= zone {
                return v % n;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn det_ln_matches_std() {
        let mut x = 1e-300;
        while x < 1e300 {
            for f in [1.0, 1.1, std::f64::consts::SQRT_2, 1.5, 1.9999, 7.3] {
                let v = x * f;
                let (a, b) = (det_ln(v), v.ln());
                assert!(
                    (a - b).abs() <= 1e-15 * b.abs().max(1.0),
                    "ln({v}): {a} vs {b}"
                );
            }
            x *= 1e7;
        }
        assert_eq!(det_ln(1.0), 0.0);
        let tiny = 5e-324;
        assert!((det_ln(tiny) - tiny.ln()).abs() < 1e-12);
    }

    #[test]
    fn streams_are_reproducible_and_distinct() {
        let draw = |seed, stream| {
            let mut r = SsaRng::new(seed, stream);
            (0..8).map(|_| r.uniform()).collect::<Vec<_>>()
        };
        assert_eq!(draw(7, 3), draw(7, 3));
        assert_ne!(draw(7, 3), draw(7, 4));
        assert_ne!(draw(7, 3), draw(8, 3));
    }

    #[test]
    fn uniform_ranges_and_moments() {
        let mut r = SsaRng::new(1, 0);
        let n = 200_000;
        let (mut sum, mut sum_e) = (0.0, 0.0);
        for _ in 0..n {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
            let v = r.uniform_open();
            assert!(v > 0.0 && v <= 1.0);
            sum += u;
            sum_e += r.exp1();
        }
        // Standard errors: 1/sqrt(12 n) = 6.5e-4 and 1/sqrt(n) = 2.2e-3.
        assert!((sum / n as f64 - 0.5).abs() < 4.0 * 6.5e-4);
        assert!((sum_e / n as f64 - 1.0).abs() < 4.0 * 2.3e-3);
    }

    #[test]
    fn below_is_uniform() {
        let mut r = SsaRng::new(2, 0);
        let mut hits = [0u32; 8];
        for _ in 0..80_000 {
            hits[r.below(8) as usize] += 1;
        }
        for h in hits {
            // Binomial(80000, 1/8): mean 10000, sd 93.5.
            assert!((h as f64 - 10_000.0).abs() < 5.0 * 93.5, "{hits:?}");
        }
    }
}
