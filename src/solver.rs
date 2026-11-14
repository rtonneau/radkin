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
    ws: &'a mut Workspace,
    sparse: Option<&'a SparseJacobian>,
}

impl OdeSystem for ChemSystem<'_> {
    fn rhs(&mut self, t: f64, y: &[f64], ydot: &mut [f64]) -> std::result::Result<(), LsodaError> {
        #[cfg(feature = "timing")]
        let started = std::time::Instant::now();
        derive(y, self.source.at(t), self.mat, self.ws, ydot);
        #[cfg(feature = "timing")]
        timing::add_rhs(started.elapsed());
        Ok(())
    }

    fn jac(&mut self, _t: f64, y: &[f64], pd: &mut [f64]) -> std::result::Result<(), LsodaError> {
        #[cfg(feature = "timing")]
        let started = std::time::Instant::now();
        jacobian(y, self.mat, self.ws, pd);
        #[cfg(feature = "timing")]
        timing::add_jac(started.elapsed());
        Ok(())
    }

    fn has_jac(&self) -> bool {
        true
    }

    fn sparsity_pattern(&self) -> Option<&SparsityPattern> {
        self.sparse.map(SparseJacobian::pattern)
    }

    fn jac_sparse(
        &mut self,
        _t: f64,
        y: &[f64],
        values: &mut [f64],
    ) -> std::result::Result<(), LsodaError> {
        let sj = self
            .sparse
            .expect("sparse Jacobian only requested in sparse mode");
        #[cfg(feature = "timing")]
        let started = std::time::Instant::now();
        jacobian_sparse(y, self.mat, sj, self.ws, values);
        #[cfg(feature = "timing")]
        timing::add_jac(started.elapsed());
        Ok(())
    }

    fn has_jac_sparse(&self) -> bool {
        self.sparse.is_some()
    }
}

/// Time breakdown of the solves run on this thread (cargo feature `timing`).
///
/// Every [`solve`] adds its kernel calls (right-hand side, dense or sparse
/// Jacobian) and its rsoda linear-solver calls ([`rsoda::TimingStats`]) to a
/// thread-local total, which [`timing::take`] returns and resets. A chunked
/// run, made of many solves, is summed.
#[cfg(feature = "timing")]
pub mod timing {
    use std::cell::Cell;
    use std::time::Duration;

    /// Call counts and wall times accumulated by [`super::solve`].
    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct Breakdown {
        /// Kernel right-hand-side calls ([`crate::kernel::derive`]).
        pub n_rhs: usize,
        pub t_rhs: Duration,
        /// Kernel Jacobian calls (dense or sparse).
        pub n_jac: usize,
        pub t_jac: Duration,
        /// rsoda `LinearSolver::prepare` calls: Jacobian evaluation (so
        /// `t_jac`, or the finite-difference RHS calls), then assembly and
        /// factorization of the iteration matrix.
        pub n_prepare: usize,
        pub t_prepare: Duration,
        /// rsoda `LinearSolver::solve` calls.
        pub n_solve: usize,
        pub t_solve: Duration,
    }

    thread_local! {
        static TOTAL: Cell<Breakdown> = Cell::new(Breakdown::default());
    }

    fn update(f: impl FnOnce(&mut Breakdown)) {
        TOTAL.with(|c| {
            let mut b = c.get();
            f(&mut b);
            c.set(b);
        });
    }

    pub(super) fn add_rhs(d: Duration) {
        update(|b| {
            b.n_rhs += 1;
            b.t_rhs += d;
        });
    }

    pub(super) fn add_jac(d: Duration) {
        update(|b| {
            b.n_jac += 1;
            b.t_jac += d;
        });
    }

    pub(super) fn add_solver(s: &rsoda::TimingStats) {
        update(|b| {
            b.n_prepare += s.n_prepare;
            b.t_prepare += s.t_prepare;
            b.n_solve += s.n_solve;
            b.t_solve += s.t_solve;
        });
    }

    /// The totals since the last call on this thread; resets them.
    pub fn take() -> Breakdown {
        TOTAL.with(|c| c.replace(Breakdown::default()))
    }
}

fn solver_failed(t_span: (f64, f64), t_reached: f64, e: &dyn std::fmt::Display) -> RadkinError {
    RadkinError::SolverFailed {
        start: t_span.0,
        stop: t_span.1,
        t_reached,
        message: e.to_string(),
    }
}

/// An rsoda solver for `system` from `(t_span.0, y0)` with the options of
/// `opts` and the first step `first_step`.
fn new_lsoda<'a>(
    system: ChemSystem<'a>,
    t_span: (f64, f64),
    y0: &[f64],
    opts: &SolveOptions,
    first_step: f64,
) -> Result<Lsoda<'static, ChemSystem<'a>>> {
    let mut builder = OptionsBuilder::new(y0.len())
        .tolerance(Tolerance::scalar(opts.rtol, opts.atol))
        .jacobian(match (opts.jacobian, opts.use_jac) {
            (JacobianMode::Dense, true) => JacobianKind::DenseUser,
            (JacobianMode::Dense, false) => JacobianKind::DenseFd,
            (JacobianMode::Sparse, true) => JacobianKind::SparseUser,
            (JacobianMode::Sparse, false) => JacobianKind::SparseFd,
        })
        .h0(first_step);
    if let Some(h) = opts.max_step {
        builder = builder.hmax(h);
    }
    let options = builder
        .build(&system)
        .map_err(|e| solver_failed(t_span, t_span.0, &e))?;
    Lsoda::new(system, options, t_span.0, y0).map_err(|e| solver_failed(t_span, t_span.0, &e))
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
    let mut ws = Workspace::new(mat);
    let sparse = (opts.jacobian == JacobianMode::Sparse).then(|| SparseJacobian::new(mat));
    let system = ChemSystem {
        mat,
        source,
        ws: &mut ws,
        sparse: sparse.as_ref(),
    };
    let mut solver = new_lsoda(system, t_span, y0, opts, opts.first_step)?;
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
            .map_err(|e| solver_failed(t_span, solver.t(), &e))?;
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
    #[cfg(feature = "timing")]
    timing::add_solver(&solver.timing());
    sol.nfev = solver.nfe();
    sol.njev = solver.nje();
    sol.nst = solver.nst();
    Ok(sol)
}

/// Kernel buffers reused by successive [`integrate_to`] calls (one per
/// thread): the [`Workspace`] and, in sparse mode, the sparsity pattern.
#[derive(Debug, Clone)]
pub struct SolveScratch {
    ws: Workspace,
    sparse: Option<SparseJacobian>,
}

impl SolveScratch {
    pub fn new(mat: &SimMatrices, mode: JacobianMode) -> Self {
        Self {
            ws: Workspace::new(mat),
            sparse: (mode == JacobianMode::Sparse).then(|| SparseJacobian::new(mat)),
        }
    }
}

/// Counters of one [`integrate_to`] and the step size to carry to the next.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Integration {
    pub nfev: usize,
    pub njev: usize,
    pub nst: usize,
    /// Last step size [s] not shortened to land on the end of the span (the
    /// last step at all if every step was).
    pub h_last: f64,
}

/// Integrate `y` [µmol/L] in place over `t_span` with the first step
/// `first_step` (`opts.first_step` and `opts.t_eval` are ignored).
///
/// Same steps as [`solve`] with `t_eval = [t_span.1]` and that first step,
/// so the result is bit-identical to its single output, without storing any
/// output and with the kernel buffers of `scratch`.
pub fn integrate_to(
    mat: &SimMatrices,
    source: &DoseRateSource,
    t_span: (f64, f64),
    y: &mut [f64],
    opts: &SolveOptions,
    first_step: f64,
    scratch: &mut SolveScratch,
) -> Result<Integration> {
    let (t0, t_end) = t_span;
    assert_eq!(
        scratch.sparse.is_some(),
        opts.jacobian == JacobianMode::Sparse,
        "scratch built for another Jacobian mode"
    );
    let system = ChemSystem {
        mat,
        source,
        ws: &mut scratch.ws,
        sparse: scratch.sparse.as_ref(),
    };
    let mut solver = new_lsoda(system, t_span, y, opts, first_step)?;
    let mut h_last = 0.0;
    let mut t = t0;
    while t < t_end {
        t = solver
            .step_tcrit(t_end, t_end)
            .map_err(|e| solver_failed(t_span, solver.t(), &e))?;
        if t < t_end || h_last == 0.0 {
            h_last = solver.hu();
        }
        if t >= t_end {
            lsoda_dense_output(&solver.nordsieck(), t, t_end, y);
        }
    }
    #[cfg(feature = "timing")]
    timing::add_solver(&solver.timing());
    Ok(Integration {
        nfev: solver.nfe(),
        njev: solver.nje(),
        nst: solver.nst(),
        h_last,
    })
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
    fn integrate_to_matches_solve_with_a_reused_scratch() {
        let cfg = Config::load(crate::consts::data_dir().join("CONFIG_PAPER_01.toml")).unwrap();
        let mut env = SimEnv::from_config(&cfg).unwrap();
        for (s, c) in [("e_aq", 2e-6), ("OH_r", 2e-6), ("H3O_plus", 2e-6)] {
            env.set_initial_cc(s, c).unwrap();
        }
        let mat = SimMatrices::build(&env).unwrap();
        let y0 = env.initial_values();
        for mode in [JacobianMode::Dense, JacobianMode::Sparse] {
            let o = SolveOptions {
                jacobian: mode,
                t_eval: Some(vec![1e-6]),
                ..opts()
            };
            let mut scratch = SolveScratch::new(&mat, mode);
            // The same scratch for several states, in both orders.
            for scale in [1.0, 3.0, 1.0] {
                let y: Vec<f64> = y0.iter().map(|v| v * scale).collect();
                let reference = solve(&mat, &DoseRateSource::None, (0.0, 1e-6), &y, &o).unwrap();
                let mut lean = y.clone();
                let r = integrate_to(
                    &mat,
                    &DoseRateSource::None,
                    (0.0, 1e-6),
                    &mut lean,
                    &o,
                    o.first_step,
                    &mut scratch,
                )
                .unwrap();
                assert_eq!(lean, reference.y[0], "{mode:?} x{scale}");
                assert_eq!(
                    (r.nfev, r.njev, r.nst),
                    (reference.nfev, reference.njev, reference.nst)
                );
                assert!(r.h_last > 0.0 && r.h_last < 1e-6);
            }
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
