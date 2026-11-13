//! Power-series polynomial with numpy's evaluation order.

/// `c[0] + c[1] x + ...` (numpy.polynomial.Polynomial).
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    pub coef: Vec<f64>,
}

impl Polynomial {
    /// Polynomial from coefficients; a single coefficient with `log10`
    /// stores `log10(c0)` (radiopyo `polynomial_factory`).
    pub fn from_coefs(coef: Vec<f64>, log10: bool) -> Self {
        if coef.len() == 1 && log10 {
            return Self {
                coef: vec![coef[0].log10()],
            };
        }
        Self { coef }
    }

    /// Constant polynomial from a scalar; with `log10` it stores `log10(v)`
    /// (`-inf` for `v <= 0`, i.e. k = 0).
    pub fn from_value(value: f64, log10: bool) -> Self {
        let c = match (log10, value >= 0.0) {
            (false, _) => value,
            (true, true) => value.log10(),
            (true, false) => f64::NEG_INFINITY,
        };
        Self { coef: vec![c] }
    }

    /// Evaluate at `x` (numpy Horner order).
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.coef.len();
        if n == 0 {
            return 0.0;
        }
        let mut c0 = self.coef[n - 1];
        for i in 2..=n {
            c0 = self.coef[n - i] + c0 * x;
        }
        c0
    }

    /// Coefficient-wise difference (zero padded).
    pub fn sub(&self, other: &Self) -> Self {
        self.zip(other, |a, b| a - b)
    }

    /// Coefficient-wise sum (zero padded).
    pub fn add(&self, other: &Self) -> Self {
        self.zip(other, |a, b| a + b)
    }

    fn zip(&self, other: &Self, op: impl Fn(f64, f64) -> f64) -> Self {
        let n = self.coef.len().max(other.coef.len());
        let at = |v: &[f64], i: usize| v.get(i).copied().unwrap_or(0.0);
        Self {
            coef: (0..n)
                .map(|i| op(at(&self.coef, i), at(&other.coef, i)))
                .collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.coef.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coef.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_uses_power_series() {
        assert_eq!(
            Polynomial::from_coefs(vec![1.0, 2.0, 3.0], false).eval(2.0),
            17.0
        );
    }

    #[test]
    fn log10_scalar_round_trips() {
        let p = Polynomial::from_value(1e10, true);
        assert_eq!(p.eval(0.003), 10.0);
        assert_eq!(10f64.powf(Polynomial::from_value(0.0, true).eval(1.0)), 0.0);
        assert_eq!(
            10f64.powf(Polynomial::from_value(-1.0, true).eval(1.0)),
            0.0
        );
    }

    #[test]
    fn subtraction_pads_coefficients() {
        let a = Polynomial::from_coefs(vec![1.0, 2.0], false);
        let b = Polynomial::from_coefs(vec![0.5], false);
        assert_eq!(a.sub(&b).coef, vec![0.5, 2.0]);
        assert_eq!(b.add(&a).coef, vec![1.5, 2.0]);
    }
}
