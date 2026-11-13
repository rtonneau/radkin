//! scipy's LSODA dense output (`scipy.integrate._ivp.lsoda`), ported exactly.
//!
//! scipy does not interpolate `t_eval` points with DINTDY. After each step,
//! `LSODA._dense_output_impl` copies the first `NQU + 1` columns of the
//! Nordsieck array (`NQU`: the order of the step just taken), rescales the
//! last one by `(HCUR / HU)^NQU` when the order is lowered for the next step
//! (`NQCUR < NQU`: ODEPACK leaves that column scaled to `HU`), and
//! `LsodaDenseOutput` evaluates the power series
//!
//! `y(t) = sum_{j=0}^{NQU} yh[:, j] * x^j`, with `x = (t - t_step) / HCUR`,
//!
//! where `t_step` is the time the step returned (`tcrit` exactly when the
//! step reached it), `x^j` is numpy's `**` (C `pow`) and the sum is
//! `np.dot`, taken here in column order.

use rsoda::Nordsieck;

/// Writes scipy's dense output at `t` into `out[..n]`, from the history `nord`
/// of the step that returned `t_step`. Allocates nothing.
///
/// # Panics
///
/// If `out` has fewer than `nord.n` components or `nord.yh` fewer than
/// `nord.nqu + 1` columns.
pub fn lsoda_dense_output(nord: &Nordsieck<'_>, t_step: f64, t: f64, out: &mut [f64]) {
    let Nordsieck {
        yh,
        nyh,
        n,
        nq,
        nqu,
        h,
        hu,
        ..
    } = *nord;
    let out = &mut out[..n];
    let col = |j: usize| &yh[j * nyh..j * nyh + n];
    let x = (t - t_step) / h;
    let rescale = (nq < nqu).then(|| (h / hu).powf(nqu as f64));
    out.copy_from_slice(col(0));
    // x ** 0 == 1 exactly: column 0 enters the dot product unchanged.
    for j in 1..=nqu {
        let p = x.powf(j as f64);
        match rescale {
            Some(r) if j == nqu => {
                for (o, &y) in out.iter_mut().zip(col(j)) {
                    *o += (y * r) * p;
                }
            }
            _ => {
                for (o, &y) in out.iter_mut().zip(col(j)) {
                    *o += y * p;
                }
            }
        }
    }
}
