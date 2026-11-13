//! LET-dependent G-values from Boscolo et al. 2020 (port of radiopyo
//! `utils/g_values.py`): a quadratic interpolating B-spline per species,
//! built as scipy `interp1d(kind="quadratic")` (= `make_interp_spline(k=2)`).

use std::path::Path;

use indexmap::IndexMap;

use crate::consts::{GF_A, GF_B};
use crate::error::{RadkinError, Result};

/// File name of the bundled Boscolo table.
pub const BOSCOLO_FILE: &str = "gvalues_from_Boscolo_2020.csv";

/// Quadratic interpolating B-spline (scipy `make_interp_spline(x, y, k=2)`).
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticSpline {
    t: Vec<f64>,
    c: Vec<f64>,
    x_min: f64,
    x_max: f64,
}

const K: usize = 2;

/// Values of the `K + 1` non-zero B-splines of degree `K` at `x` in knot
/// interval `l` (`t[l] <= x < t[l+1]`), Cox–de Boor recursion.
fn basis(t: &[f64], l: usize, x: f64) -> [f64; K + 1] {
    let mut n = [0.0; K + 1];
    n[0] = 1.0;
    let mut left = [0.0; K + 1];
    let mut right = [0.0; K + 1];
    for j in 1..=K {
        left[j] = x - t[l + 1 - j];
        right[j] = t[l + j] - x;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

impl QuadraticSpline {
    /// Interpolant of `(x, y)`, `x` strictly increasing, at least 3 points.
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        let n = x.len();
        if n < K + 1 || y.len() != n {
            return Err(RadkinError::GValue(
                "spline needs >= 3 points and len(x) == len(y)".into(),
            ));
        }
        // Knots: Greville-like sites, skipping the 2nd and 2nd-to-last ones.
        let mut t = vec![x[0]; K + 1];
        let mids: Vec<f64> = (0..n - 1).map(|i| (x[i + 1] + x[i]) / 2.0).collect();
        t.extend_from_slice(&mids[1..mids.len() - 1]);
        t.extend(std::iter::repeat_n(x[n - 1], K + 1));
        // Collocation matrix, dense LU with partial pivoting.
        let mut a = vec![0.0; n * n];
        for (i, &xi) in x.iter().enumerate() {
            let l = Self::interval(&t, n, xi);
            let b = basis(&t, l, xi);
            for (r, v) in b.iter().enumerate() {
                a[i * n + (l - K + r)] = *v;
            }
        }
        let c = solve(&mut a, y.to_vec(), n)?;
        Ok(Self {
            t,
            c,
            x_min: x[0],
            x_max: x[n - 1],
        })
    }

    fn interval(t: &[f64], n: usize, x: f64) -> usize {
        // Largest l in [K, n-1] with t[l] <= x (x == t[n] uses l = n - 1).
        let mut l = K;
        while l + 1 < n && t[l + 1] <= x {
            l += 1;
        }
        l
    }

    /// Value at `x`; `Err` outside `[x0, xN]` (interp1d `bounds_error=True`).
    pub fn eval(&self, x: f64) -> Result<f64> {
        if !(x >= self.x_min && x <= self.x_max) {
            return Err(RadkinError::GValue(format!(
                "LET {x} keV/µm outside the Boscolo data range [{}, {}]",
                self.x_min, self.x_max
            )));
        }
        let n = self.c.len();
        let l = Self::interval(&self.t, n, x);
        let b = basis(&self.t, l, x);
        Ok((0..=K).map(|r| self.c[l - K + r] * b[r]).sum())
    }
}

fn solve(a: &mut [f64], mut b: Vec<f64>, n: usize) -> Result<Vec<f64>> {
    for col in 0..n {
        let p = (col..n)
            .max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))
            .unwrap();
        if a[p * n + col] == 0.0 {
            return Err(RadkinError::GValue("singular collocation matrix".into()));
        }
        if p != col {
            for k in 0..n {
                a.swap(p * n + k, col * n + k);
            }
            b.swap(p, col);
        }
        for i in col + 1..n {
            let f = a[i * n + col] / a[col * n + col];
            if f != 0.0 {
                for k in col..n {
                    a[i * n + k] -= f * a[col * n + k];
                }
                b[i] -= f * b[col];
            }
        }
    }
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|k| a[i * n + k] * b[k]).sum();
        b[i] = (b[i] - s) / a[i * n + i];
    }
    Ok(b)
}

/// Boscolo 2020 G-values [#/100 eV] versus LET [keV/µm].
#[derive(Debug, Clone, PartialEq)]
pub struct GValuesBoscolo {
    splines: IndexMap<String, QuadraticSpline>,
}

impl GValuesBoscolo {
    /// Load `gvalues_from_Boscolo_2020.csv` from `data_dir`.
    pub fn load(data_dir: &Path) -> Result<Self> {
        let path = data_dir.join(BOSCOLO_FILE);
        let text = std::fs::read_to_string(&path).map_err(|e| RadkinError::io(&path, e))?;
        Self::from_csv_str(&text, &path.display().to_string())
    }

    /// Parse the text of a Boscolo CSV (`source` names it in errors).
    pub fn from_csv_str(text: &str, source: &str) -> Result<Self> {
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header: Vec<String> = lines
            .next()
            .unwrap_or_default()
            .split(',')
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        let mut x = Vec::new();
        let mut cols: Vec<Vec<f64>> = vec![Vec::new(); header.len()];
        for line in lines {
            let vals: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>())
                .collect::<std::result::Result<_, _>>()
                .map_err(|e| RadkinError::GValue(format!("{source}: {e}")))?;
            x.push(vals[0]);
            for (c, v) in cols.iter_mut().zip(&vals[1..]) {
                c.push(*v);
            }
        }
        let splines = header
            .into_iter()
            .zip(cols)
            .map(|(h, y)| Ok((h, QuadraticSpline::new(&x, &y)?)))
            .collect::<Result<_>>()?;
        Ok(Self { splines })
    }

    /// Whether `species` has Boscolo data.
    pub fn contains(&self, species: &str) -> bool {
        self.splines.contains_key(species)
    }

    /// Dose-rate factor `0.97207704 * dr^-0.10972` (1 if `None`).
    pub fn scale_factor(dose_rate: Option<f64>) -> f64 {
        match dose_rate {
            None => 1.0,
            Some(dr) => GF_A * dr.powf(GF_B),
        }
    }

    /// G-value of `species` at `let_`; H3O_plus takes the e_aq value (ADR 0002).
    pub fn value_at(&self, species: &str, let_: f64, dose_rate: Option<f64>) -> Result<f64> {
        let species = if species == "H3O_plus" {
            "e_aq"
        } else {
            species
        };
        let s = self
            .splines
            .get(species)
            .ok_or_else(|| RadkinError::GValue(format!("no Boscolo data for '{species}'")))?;
        Ok(s.eval(let_)? * Self::scale_factor(dose_rate))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gv() -> GValuesBoscolo {
        GValuesBoscolo::load(&crate::consts::data_dir()).unwrap()
    }

    #[test]
    fn spline_interpolates_the_data_points() {
        let g = gv();
        assert!((g.value_at("e_aq", 0.13, None).unwrap() - 2.7).abs() < 1e-12);
        assert!((g.value_at("OH_r", 0.5, None).unwrap() - 2.765496688478056).abs() < 1e-12);
    }

    #[test]
    fn h3o_plus_equals_e_aq() {
        let g = gv();
        assert_eq!(
            g.value_at("H3O_plus", 3.2, None).unwrap(),
            g.value_at("e_aq", 3.2, None).unwrap()
        );
    }

    #[test]
    fn let_out_of_range_is_an_error() {
        assert!(gv().value_at("OH_r", 200.0, None).is_err());
        assert!(gv().value_at("OH_r", 0.1, None).is_err());
    }
}
