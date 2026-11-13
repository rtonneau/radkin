//! `scipy.integrate.solve_ivp(method="LSODA")` equivalent on rsoda.
//!
//! scipy drives LSODA one step at a time (`ITASK = 5`, `TCRIT = t_bound`)
//! and fills `t_eval` from the Nordsieck history of each step with its own
//! dense output (not DINTDY). Here: [`rsoda::Lsoda::step_tcrit`] for the
//! steps, and [`crate::dense::lsoda_dense_output`] on
//! [`rsoda::Lsoda::nordsieck`] for the `t_eval` points of each step,
//! including one equal to the step's end, as scipy does.

use rsoda::{
    JacobianKind, Lsoda, LsodaError, OdeSystem, OptionsBuilder, SparsityPattern, Tolerance,
};

use crate::dense::lsoda_dense_output;
use crate::error::{RadkinError, Result};
use crate::kernel::{SparseJacobian, Workspace, derive, jacobian, jacobian_sparse};
use crate::matrices::SimMatrices;
use crate::physics::Beam;

/// Dose rate seen by the right-hand side.
#[derive(Debug, Clone, PartialEq)]
pub enum DoseRateSource {
    /// No irradiation.
    None,
    /// Constant dose rate [Gy/s] (one chunk of a chunked run).
    Constant(f64),
    /// The beam, read at every evaluation.
    Beam(Beam),
}

impl DoseRateSource {
    pub fn at(&self, t: f64) -> f64 {
        match self {
            Self::None => 0.0,
            Self::Constant(d) => *d,
            Self::Beam(b) => b.at(t),
        }
    }
}

/// Storage of the Jacobian and linear solver of the stiff (BDF) steps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum JacobianMode {
    /// Dense LU (rsoda's LINPACK port, as scipy's LSODA).
    #[default]
    Dense,
    /// Sparse LU on the kernel's sparsity pattern
    /// ([`crate::kernel::SparseJacobian`]).
    Sparse,
}

/// Step options (scipy `solve_ivp` keyword arguments).
#[derive(Debug, Clone, PartialEq)]
pub struct SolveOptions {
    pub rtol: f64,
    pub atol: f64,
    pub first_step: f64,
    pub max_step: Option<f64>,
    /// Analytic Jacobian (else finite differences).
    pub use_jac: bool,
    /// Dense or sparse Jacobian and linear solver.
    pub jacobian: JacobianMode,
    /// Output times; `None` stores every step (and the initial point).
    pub t_eval: Option<Vec<f64>>,
}

/// Solution of one span.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Solution {
    pub t: Vec<f64>,
    /// State at each time [µmol/L].
    pub y: Vec<Vec<f64>>,
    pub nfev: usize,
    pub njev: usize,
    pub nst: usize,
}

struct ChemSystem<'a> {
    mat: &'a SimMatrices,
    source: &'a DoseRateSource,
    ws: Workspace,
    sparse: Option<SparseJacobian>,
}

impl OdeSystem for ChemSystem<'_> {
    fn rhs(&mut self, t: f64, y: &[f64], ydot: &mut [f64]) -> std::result::Result<(), LsodaError> {
        derive(y, self.source.at(t), self.mat, &mut self.ws, ydot);
        Ok(())
    }

    fn jac(&mut self, _t: f64, y: &[f64], pd: &mut [f64]) -> std::result::Result<(), LsodaError> {
        jacobian(y, self.mat, &mut self.ws, pd);
        Ok(())
    }

    fn has_jac(&self) -> bool {
        true
    }

    fn sparsity_pattern(&self) -> Option<&SparsityPattern> {
        self.sparse.as_ref().map(SparseJacobian::pattern)
    }

    fn jac_sparse(
        &mut self,
        _t: f64,
        y: &[f64],
        values: &mut [f64],
    ) -> std::result::Result<(), LsodaError> {
        let sj = self
            .sparse
            .as_ref()
            .expect("sparse Jacobian only requested in sparse mode");
        jacobian_sparse(y, self.mat, sj, &mut self.ws, values);
        Ok(())
    }

    fn has_jac_sparse(&self) -> bool {
        self.sparse.is_some()
    }
}

/// Integrate the dynamic state `y0` [µmol/L] over `t_span`.
pub fn solve(
    mat: &SimMatrices,
    source: &DoseRateSource,
    t_span: (f64, f64),
    y0: &[f64],
    opts: &SolveOptions,
) -> Result<Solution> {
    let (t0, t_end) = t_span;
    let fail = |t_reached: f64, e: &dyn std::fmt::Display| RadkinError::SolverFailed {
        start: t0,
        stop: t_end,
        t_reached,
        message: e.to_string(),
    };
    if let Some(te) = &opts.t_eval {
        if te.iter().any(|&t| t < t0 || t > t_end) {
            return Err(RadkinError::Value(format!(
                "Values in `t_eval` are not within `t_span` ({t0}, {t_end})"
            )));
        }
        if te.windows(2).any(|w| w[1] < w[0]) {
            return Err(RadkinError::Value(
                "Values in `t_eval` are not properly sorted".into(),
            ));
        }
    }
    let system = ChemSystem {
        mat,
        source,
        ws: Workspace::new(mat),
        sparse: (opts.jacobian == JacobianMode::Sparse).then(|| SparseJacobian::new(mat)),
    };
    let mut builder = OptionsBuilder::new(y0.len())
        .tolerance(Tolerance::scalar(opts.rtol, opts.atol))
        .jacobian(match (opts.jacobian, opts.use_jac) {
            (JacobianMode::Dense, true) => JacobianKind::DenseUser,
            (JacobianMode::Dense, false) => JacobianKind::DenseFd,
            (JacobianMode::Sparse, true) => JacobianKind::SparseUser,
            (JacobianMode::Sparse, false) => JacobianKind::SparseFd,
        })
        .h0(opts.first_step);
    if let Some(h) = opts.max_step {
        builder = builder.hmax(h);
    }
    let options = builder.build(&system).map_err(|e| fail(t0, &e))?;
    let mut solver = Lsoda::new(system, options, t0, y0).map_err(|e| fail(t0, &e))?;
    let mut sol = Solution::default();
    let t_eval = opts.t_eval.as_deref();
    if t_eval.is_none() {
        sol.t.push(t0);
        sol.y.push(y0.to_vec());
    }
    let mut next_eval = 0;
    let mut t = t0;
    while t < t_end {
        t = solver
            .step_tcrit(t_end, t_end)
            .map_err(|e| fail(solver.t(), &e))?;
        match t_eval {
            None => {
                sol.t.push(t);
                sol.y.push(solver.y().to_vec());
            }
            Some(te) => {
                let nord = solver.nordsieck();
                while next_eval < te.len() && te[next_eval] <= t {
                    let tout = te[next_eval];
                    let mut y = vec![0.0; y0.len()];
                    lsoda_dense_output(&nord, t, tout, &mut y);
                    sol.t.push(tout);
                    sol.y.push(y);
                    next_eval += 1;
                }
            }
        }
    }
    sol.nfev = solver.nfe();
    sol.njev = solver.nje();
    sol.nst = solver.nst();
    Ok(sol)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::matrices::SimMatrices;

    fn decay() -> (SimMatrices, Vec<f64>, usize) {
        let cfg = Config::from_toml_str(
            "[parameters]\n[concentrations.initial]\nA = 1e-3\n[[reactions.k_reaction]]\nreaction = \"A -> B\"\nk_value = 1.0\n",
            None,
        )
        .unwrap();
        let env = SimEnv::from_config(&cfg).unwrap();
        let ia = env.species.index_of("A").unwrap();
        (SimMatrices::build(&env).unwrap(), env.initial_values(), ia)
    }

    fn opts() -> SolveOptions {
        SolveOptions {
            rtol: 1e-10,
            atol: 1e-12,
            first_step: 1e-15,
            max_step: None,
            use_jac: true,
            jacobian: JacobianMode::Dense,
            t_eval: None,
        }
    }

    #[test]
    fn exponential_decay_reaches_the_analytic_value() {
        let (mat, y0, ia) = decay();
        let sol = solve(&mat, &DoseRateSource::None, (0.0, 1.0), &y0, &opts()).unwrap();
        assert_eq!(sol.t[0], 0.0);
        assert_eq!(*sol.t.last().unwrap(), 1.0);
        let a = sol.y.last().unwrap()[ia];
        assert!((a - 1000.0 * (-1.0f64).exp()).abs() < 1e-8 * 1000.0, "{a}");
        assert!(sol.nfev > 0);
    }

    #[test]
    fn t_eval_points_are_interpolated() {
        let (mat, y0, ia) = decay();
        let o = SolveOptions {
            t_eval: Some(vec![0.0, 0.5, 1.0]),
            use_jac: false,
            jacobian: JacobianMode::Dense,
            ..opts()
        };
        let sol = solve(&mat, &DoseRateSource::None, (0.0, 1.0), &y0, &o).unwrap();
        assert_eq!(sol.t, vec![0.0, 0.5, 1.0]);
        for (t, y) in sol.t.iter().zip(&sol.y) {
            assert!(
                (y[ia] - 1000.0 * (-t).exp()).abs() < 1e-6,
                "t={t}: {}",
                y[ia]
            );
        }
    }

    #[test]
    fn t_eval_outside_span_is_an_error() {
        let (mat, y0, _) = decay();
        let o = SolveOptions {
            t_eval: Some(vec![0.5, 2.0]),
            ..opts()
        };
        assert!(solve(&mat, &DoseRateSource::None, (0.0, 1.0), &y0, &o).is_err());
    }

    #[test]
    fn solver_failure_is_an_error_not_a_partial_result() {
        let (mat, y0, _) = decay();
        let o = SolveOptions {
            rtol: 0.0,
            atol: 0.0,
            ..opts()
        };
        let err = solve(&mat, &DoseRateSource::None, (0.0, 1.0), &y0, &o).unwrap_err();
        assert!(
            matches!(err, crate::RadkinError::SolverFailed { .. }),
            "{err}"
        );
    }
}
