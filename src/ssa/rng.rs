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

/// `exp(x)` from IEEE basic operations only (relative error below 4e-16).
///
/// `x = k ln 2 + r` with `|r| <= ln(2) / 2`, then the Taylor series of
/// `exp(r)` and a scaling by `2^k`.
// The published digits are kept, beyond what an f64 holds.
#[allow(clippy::excessive_precision)]
pub fn det_exp(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > 709.782_712_893_384 {
        return f64::INFINITY;
    }
    if x < -745.2 {
        return 0.0;
    }
    /// High and low parts of ln 2 (fdlibm).
    const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-01;
    const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
    let k = (x * std::f64::consts::LOG2_E).round();
    let r = (x - k * LN2_HI) - k * LN2_LO;
    // 1 + r + r^2/2! + ... + r^17/17! (Horner); the next term is below 3e-24.
    let mut p = 1.0;
    let mut n = 17.0;
    while n >= 1.0 {
        p = 1.0 + p * r / n;
        n -= 1.0;
    }
    // 2^k in two factors, so that a subnormal result is rounded once.
    let k = k as i64;
    let half = k / 2;
    let pow2 = |e: i64| f64::from_bits(((e + 1023) as u64) << 52);
    p * pow2(half) * pow2(k - half)
}

/// Quantile of the standard normal distribution at `p` in `(0, 1)`
/// (Wichura, algorithm AS 241, PPND16; relative error about 1e-16), from
/// IEEE basic operations and [`det_ln`]. `G4ErrorFunction::NormQuantile` of
/// Geant4 is the same algorithm.
// The published digits are kept, beyond what an f64 holds.
#[allow(clippy::excessive_precision)]
pub fn norm_quantile(p: f64) -> f64 {
    const A: [f64; 8] = [
        3.387_132_872_796_366_608,
        1.331_416_678_917_843_774_5e2,
        1.971_590_950_306_551_442_7e3,
        1.373_169_376_550_946_112_5e4,
        4.592_195_393_154_987_145_7e4,
        6.726_577_092_700_870_085_3e4,
        3.343_057_558_358_812_810_5e4,
        2.509_080_928_730_122_672_7e3,
    ];
    const B: [f64; 8] = [
        1.0,
        4.231_333_070_160_091_125_2e1,
        6.871_870_074_920_579_083e2,
        5.394_196_021_424_751_107_7e3,
        2.121_379_430_158_659_586_7e4,
        3.930_789_580_009_271_061e4,
        2.872_908_573_572_194_267_4e4,
        5.226_495_278_852_854_561e3,
    ];
    const C: [f64; 8] = [
        1.423_437_110_749_683_577_34,
        4.630_337_846_156_545_295_9,
        5.769_497_221_460_691_405_5,
        3.647_848_324_763_204_605_04,
        1.270_458_252_452_368_382_58,
        2.417_807_251_774_506_117_7e-1,
        2.272_384_498_926_918_458_33e-2,
        7.745_450_142_783_414_076_4e-4,
    ];
    const D: [f64; 8] = [
        1.0,
        2.053_191_626_637_758_821_87,
        1.676_384_830_183_803_849_4,
        6.897_673_349_851_000_045_5e-1,
        1.481_039_764_274_800_745_9e-1,
        1.519_866_656_361_645_719_66e-2,
        5.475_938_084_995_344_946e-4,
        1.050_750_071_644_416_843_24e-9,
    ];
    const E: [f64; 8] = [
        6.657_904_643_501_103_777_2,
        5.463_784_911_164_114_369_9,
        1.784_826_539_917_291_335_8,
        2.965_605_718_285_048_912_3e-1,
        2.653_218_952_657_612_309_3e-2,
        1.242_660_947_388_078_438_6e-3,
        2.711_555_568_743_487_578_15e-5,
        2.010_334_399_292_288_132_65e-7,
    ];
    const F: [f64; 8] = [
        1.0,
        5.998_322_065_558_879_376_9e-1,
        1.369_298_809_227_358_053_1e-1,
        1.487_536_129_085_061_485_25e-2,
        7.868_691_311_456_132_591e-4,
        1.846_318_317_510_054_681_8e-5,
        1.421_511_758_316_445_888_7e-7,
        2.044_263_103_389_939_785_64e-15,
    ];
    fn horner(c: &[f64; 8], r: f64) -> f64 {
        c.iter().rev().fold(0.0, |acc, &x| acc * r + x)
    }
    let q = p - 0.5;
    if q.abs() < 0.425 {
        let r = 0.180_625 - q * q;
        return q * horner(&A, r) / horner(&B, r);
    }
    let tail = if q < 0.0 { p } else { 1.0 - p };
    if tail <= 0.0 {
        return if q < 0.0 {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    let r = (-det_ln(tail)).sqrt();
    let z = if r <= 5.0 {
        let r = r - 1.6;
        horner(&C, r) / horner(&D, r)
    } else {
        let r = r - 5.0;
        horner(&E, r) / horner(&F, r)
    };
    if q < 0.0 { -z } else { z }
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

    /// Uniform in `(0, 1)`, both ends excluded.
    pub fn uniform_inside(&mut self) -> f64 {
        ((self.0.next_u64() >> 11) as f64 + 0.5) * ULP53
    }

    /// Standard normal, by inversion ([`norm_quantile`]): one uniform each.
    pub fn normal(&mut self) -> f64 {
        norm_quantile(self.uniform_inside())
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

    #[test]
    fn det_exp_matches_std() {
        let mut x = -708.0;
        while x < 709.7 {
            let (a, b) = (det_exp(x), x.exp());
            assert!((a - b).abs() <= 4e-16 * b, "exp({x}): {a} against {b}");
            x += 0.731;
        }
        for x in [0.0, 1e-300, -1e-300, 1e-9, -1e-9, 0.3465, -0.3465, 1.0] {
            let (a, b) = (det_exp(x), x.exp());
            assert!((a - b).abs() <= 4e-16 * b, "exp({x}): {a} against {b}");
        }
        // subnormal results
        for x in [-710.0, -730.0, -744.0] {
            let (a, b) = (det_exp(x), x.exp());
            assert!(
                (a - b).abs() <= 1e-10 * b + 5e-324,
                "exp({x}): {a} against {b}"
            );
        }
        assert_eq!(det_exp(0.0), 1.0);
        assert_eq!(det_exp(-800.0), 0.0);
        assert_eq!(det_exp(710.0), f64::INFINITY);
        assert!(det_exp(f64::NAN).is_nan());
    }

    #[test]
    fn norm_quantile_known_values() {
        // quantiles of the standard normal distribution
        let cases = [
            (0.5, 0.0),
            (0.975, 1.959_963_984_540_054_2),
            (0.001, -3.090_232_306_167_813_5),
            (1e-10, -6.361_340_902_404_056),
            (1e-100, -21.273_453_560_965_324),
        ];
        for (p, z) in cases {
            let got = norm_quantile(p);
            assert!((got - z).abs() <= 2e-15 * z.abs().max(1.0), "{p}: {got}");
        }
        assert_eq!(norm_quantile(0.0), f64::NEG_INFINITY);
        assert_eq!(norm_quantile(1.0), f64::INFINITY);
        // antisymmetric
        assert_eq!(norm_quantile(0.25), -norm_quantile(0.75));
    }

    #[test]
    fn normal_moments() {
        let mut rng = SsaRng::new(11, 0);
        let n = 400_000;
        let (mut s1, mut s2, mut s4) = (0.0, 0.0, 0.0);
        for _ in 0..n {
            let z = rng.normal();
            assert!(z.is_finite());
            s1 += z;
            s2 += z * z;
            s4 += z * z * z * z;
        }
        let n = n as f64;
        assert!((s1 / n).abs() < 4.0 / n.sqrt(), "mean {}", s1 / n);
        assert!(
            (s2 / n - 1.0).abs() < 4.0 * (2.0 / n).sqrt(),
            "var {}",
            s2 / n
        );
        assert!(
            (s4 / n - 3.0).abs() < 4.0 * (96.0 / n).sqrt(),
            "m4 {}",
            s4 / n
        );
    }
}
