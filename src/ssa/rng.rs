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

/// `ln(k!)` from [`det_ln`]: a sum of logarithms below 20, the Stirling
/// series above (absolute error below 1e-13).
pub fn ln_factorial(k: u64) -> f64 {
    if k < 20 {
        return (2..=k).map(|i| det_ln(i as f64)).sum();
    }
    /// ln(2 pi) / 2.
    const HALF_LN_2PI: f64 = 0.918_938_533_204_672_7;
    let x = k as f64;
    let x2 = x * x;
    let series = (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0 - 1.0 / (1680.0 * x2)) / x2) / x2) / x;
    (x + 0.5) * det_ln(x) - x + HALF_LN_2PI + series
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

    /// Poisson count of mean `lambda` (0 for `lambda <= 0`).
    ///
    /// Below 10: the number of unit-rate exponential arrivals before
    /// `lambda`. From 10 on: Hörmann's transformed rejection (PTRS, 1993).
    pub fn poisson(&mut self, lambda: f64) -> i64 {
        if lambda.is_nan() || lambda <= 0.0 {
            return 0;
        }
        if lambda < 10.0 {
            let mut k = 0;
            let mut arrival = self.exp1();
            while arrival < lambda {
                k += 1;
                arrival += self.exp1();
            }
            return k;
        }
        let ln_lambda = det_ln(lambda);
        let b = 0.931 + 2.53 * lambda.sqrt();
        let a = -0.059 + 0.024_83 * b;
        let inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
        let v_r = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = self.uniform() - 0.5;
            let v = self.uniform_open();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor();
            if us >= 0.07 && v <= v_r {
                return k as i64;
            }
            if k < 0.0 || (us < 0.013 && v > us) {
                continue;
            }
            let lhs = det_ln(v) + det_ln(inv_alpha) - det_ln(a / (us * us) + b);
            if lhs <= k * ln_lambda - lambda - ln_factorial(k as u64) {
                return k as i64;
            }
        }
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
    fn ln_factorial_matches_the_sum_of_logs() {
        let mut exact = 0.0f64;
        for k in 1..=400u64 {
            exact += (k as f64).ln();
            let got = ln_factorial(k);
            assert!(
                (got - exact).abs() <= 1e-12 * exact.max(1.0),
                "{k}: {got} vs {exact}"
            );
        }
        assert_eq!(ln_factorial(0), 0.0);
        assert_eq!(ln_factorial(1), 0.0);
    }

    #[test]
    fn poisson_mean_and_variance() {
        let mut r = SsaRng::new(3, 0);
        assert_eq!(r.poisson(0.0), 0);
        assert_eq!(r.poisson(-1.0), 0);
        let n = 100_000;
        for lambda in [0.05, 0.9, 4.0, 9.99, 10.0, 12.5, 37.0, 1e3, 2.5e6] {
            let (mut sum, mut sum_sq) = (0.0, 0.0);
            for _ in 0..n {
                let k = r.poisson(lambda);
                assert!(k >= 0);
                let d = k as f64 - lambda;
                sum += d;
                sum_sq += d * d;
            }
            let nf = n as f64;
            let mean = sum / nf;
            let var = sum_sq / nf - mean * mean;
            // sd of the mean: sqrt(lambda / n); of the variance: about
            // lambda sqrt((2 + 1 / lambda) / n).
            assert!(
                mean.abs() <= 4.5 * (lambda / nf).sqrt(),
                "lambda {lambda}: mean off by {mean}"
            );
            let sd_var = lambda * ((2.0 + 1.0 / lambda) / nf).sqrt();
            assert!(
                (var - lambda).abs() <= 4.5 * sd_var,
                "lambda {lambda}: variance {var}"
            );
        }
    }

    #[test]
    fn poisson_follows_the_probabilities() {
        // Chi-square of the counts against the exact probabilities, on both
        // sides of the switch at 10.
        let mut r = SsaRng::new(4, 0);
        let n = 200_000;
        for lambda in [3.2f64, 15.0, 60.0] {
            let k_max = (lambda + 8.0 * lambda.sqrt()) as usize + 2;
            let mut hits = vec![0u32; k_max + 1];
            for _ in 0..n {
                let k = (r.poisson(lambda) as usize).min(k_max);
                hits[k] += 1;
            }
            let mut ln_p = -lambda;
            let (mut chi2, mut dof, mut tail_p, mut tail_hits) = (0.0, 0usize, 1.0f64, 0u32);
            let mut pending = (0.0f64, 0u32);
            for (k, h) in hits.iter().enumerate().take(k_max) {
                if k > 0 {
                    ln_p += lambda.ln() - (k as f64).ln();
                }
                let p = ln_p.exp();
                tail_p -= p;
                pending.0 += p * n as f64;
                pending.1 += h;
                // Pool neighbours until at least 20 are expected.
                if pending.0 >= 20.0 {
                    chi2 += (pending.1 as f64 - pending.0).powi(2) / pending.0;
                    dof += 1;
                    pending = (0.0, 0);
                }
            }
            tail_hits += hits[k_max] + pending.1;
            let tail_expected = tail_p.max(0.0) * n as f64 + pending.0;
            if tail_expected >= 5.0 {
                chi2 += (tail_hits as f64 - tail_expected).powi(2) / tail_expected;
                dof += 1;
            }
            let dof = (dof - 1) as f64;
            // chi2 has mean dof and sd sqrt(2 dof).
            assert!(
                chi2 <= dof + 5.0 * (2.0 * dof).sqrt(),
                "lambda {lambda}: chi2 {chi2} for {dof} degrees of freedom"
            );
        }
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
