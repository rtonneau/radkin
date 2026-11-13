//! Per-cell ODE kernel: reaction rates, right-hand side and analytic
//! Jacobian (port of radiopyo `simulation/ode_functions_solver.py`).
//!
//! Every function works on one cell's state slice and caller-provided
//! buffers (no allocation per call), so a voxel driver can call it per cell.
//!
//! - Rate laws read the clipped concentrations `max(c, 0)` (ADR 0004).
//! - Each "negative" product multiplies the rate by `c⁺/(c⁺+ε)`,
//!   `ε = DEPLETION_LIMIT` (ADR 0005).
//! - A dynamic species below 0 is read as 0 by every rate, its own row is
//!   continued linearly from `c = 0⁺` (overshoot relaxation, ADR 0007).

use crate::consts::DEPLETION_LIMIT;
use rsoda::SparsityPattern;

use crate::matrices::{EzOvershootEntry, KOvershootEntry, SimMatrices};

const EPS: f64 = DEPLETION_LIMIT;

/// Gradient accumulator and scratch list of the overshoot path.
type GradBuffers<'a> = Option<(&'a mut Vec<f64>, &'a mut Vec<(usize, f64)>)>;

/// Rates of every reaction [mol/L/s].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReactionRates {
    pub k: Vec<f64>,
    pub radiolytic: Vec<f64>,
    pub enzymatic: Vec<f64>,
}

impl ReactionRates {
    pub fn new(mat: &SimMatrices) -> Self {
        Self {
            k: vec![0.0; mat.k.k_values.len()],
            radiolytic: vec![0.0; mat.radiolytic.yields.len()],
            enzymatic: vec![0.0; mat.enzymatic.k_values.len()],
        }
    }

    /// k, radiolytic, enzymatic (radiopyo `ReactionRates.stacked`).
    pub fn stacked(&self) -> Vec<f64> {
        let mut v = self.k.clone();
        v.extend(&self.radiolytic);
        v.extend(&self.enzymatic);
        v
    }
}

/// Per-species values derived from a state, padded with the factor-1 slot.
#[derive(Debug, Clone)]
struct Prepared {
    raw: Vec<f64>,
    clip: Vec<f64>,
    /// `c⁺/(c⁺+ε)`.
    sw: Vec<f64>,
    /// `d sw / dc` (0 below 0).
    sl: Vec<f64>,
}

impl Prepared {
    fn new(n_species: usize) -> Self {
        Self {
            raw: vec![0.0; n_species],
            clip: vec![1.0; n_species + 1],
            sw: vec![1.0; n_species + 1],
            sl: vec![0.0; n_species + 1],
        }
    }

    fn set(&mut self, i: usize, v: f64) {
        self.raw[i] = v;
        let c = v.max(0.0);
        self.clip[i] = c;
        self.sw[i] = c / (c + EPS);
        self.sl[i] = if v < 0.0 {
            0.0
        } else {
            EPS / ((c + EPS) * (c + EPS))
        };
    }

    fn fill(&mut self, y_full: &[f64]) {
        for (i, v) in y_full.iter().enumerate() {
            self.set(i, *v);
        }
    }

    fn copy_from(&mut self, o: &Prepared) {
        self.raw.copy_from_slice(&o.raw);
        self.clip.copy_from_slice(&o.clip);
        self.sw.copy_from_slice(&o.sw);
        self.sl.copy_from_slice(&o.sl);
    }
}

/// Buffers reused across calls.
#[derive(Debug, Clone)]
pub struct Workspace {
    y_full: Vec<f64>,
    p: Prepared,
    p0: Prepared,
    rates: ReactionRates,
    acc: Vec<f64>,
    deriv: Vec<(usize, f64)>,
    grad: Vec<f64>,
}

impl Workspace {
    pub fn new(mat: &SimMatrices) -> Self {
        Self {
            y_full: vec![0.0; mat.n_species],
            p: Prepared::new(mat.n_species),
            p0: Prepared::new(mat.n_species),
            rates: ReactionRates::new(mat),
            acc: vec![0.0; mat.n_dyn],
            deriv: Vec::with_capacity(16),
            grad: vec![0.0; mat.n_species],
        }
    }

    /// Rates computed by the last [`derive`] call.
    pub fn rates(&self) -> &ReactionRates {
        &self.rates
    }
}

/// Full molar state: `y` [µmol/L] (dynamic) then the fixed species.
pub fn to_full_molar(y: &[f64], mat: &SimMatrices, out: &mut [f64]) {
    let n = mat.n_dyn;
    for (o, v) in out[..n].iter_mut().zip(y) {
        *o = v * 1e-6;
    }
    out[n..].copy_from_slice(&mat.cst_sp);
}

fn mass_action(clip: &[f64], k: f64, slots: &[usize]) -> f64 {
    let mut prod = 1.0;
    for &s in slots {
        prod *= clip[s];
    }
    prod * k
}

fn switch(sw: &[f64], negative: &[usize]) -> f64 {
    let mut prod = 1.0;
    for &s in negative {
        prod *= sw[s];
    }
    prod
}

/// `d(mass_action)/dc_s` pushed into `out` (columns of negative species are 0).
fn mass_action_derivatives(
    p: &Prepared,
    k: f64,
    slots: &[usize],
    out: &mut Vec<(usize, f64)>,
    scale: f64,
) {
    let n_species = p.raw.len();
    for (a, &s) in slots.iter().enumerate() {
        if s >= n_species || p.raw[s] < 0.0 {
            continue;
        }
        let mut others = 1.0;
        for (b, &q) in slots.iter().enumerate() {
            if b != a {
                others *= p.clip[q];
            }
        }
        out.push((s, others * k * scale));
    }
}

/// `d(switch)/dc_s` pushed into `out`.
fn switch_derivatives(p: &Prepared, negative: &[usize], out: &mut Vec<(usize, f64)>, scale: f64) {
    for (a, &s) in negative.iter().enumerate() {
        let mut others = 1.0;
        for (b, &q) in negative.iter().enumerate() {
            if b != a {
                others *= p.sw[q];
            }
        }
        out.push((s, others * p.sl[s] * scale));
    }
}

fn michaelis(clip: &[f64], k: f64, km: f64, enzyme: usize, substrate: usize) -> f64 {
    let e = clip[enzyme];
    let s = clip[substrate];
    k * e * s / (km + s)
}

fn compute_rates(p: &Prepared, dose_rate: f64, mat: &SimMatrices, out: &mut ReactionRates) {
    let k = &mat.k;
    for r in 0..k.k_values.len() {
        let law = mass_action(&p.clip, k.k_values[r], &k.reactant_index[r]);
        out.k[r] = if k.negative_index[r].is_empty() {
            law
        } else {
            law * switch(&p.sw, &k.negative_index[r])
        };
    }
    for (o, y) in out.radiolytic.iter_mut().zip(&mat.radiolytic.yields) {
        *o = y * dose_rate;
    }
    let ez = &mat.enzymatic;
    for r in 0..ez.k_values.len() {
        let law = michaelis(
            &p.clip,
            ez.k_values[r],
            ez.k_michaelis[r],
            ez.enzyme_index[r],
            ez.substrate_index[r],
        );
        out.enzymatic[r] = if ez.negative_index[r].is_empty() {
            law
        } else {
            law * switch(&p.sw, &ez.negative_index[r])
        };
    }
}

/// Rates of every reaction for a full molar state `y_full` (fixed species
/// included) at `dose_rate` [Gy/s] (radiopyo `reaction_rates`).
pub fn reaction_rates(
    y_full: &[f64],
    dose_rate: f64,
    mat: &SimMatrices,
    ws: &mut Workspace,
    out: &mut ReactionRates,
) {
    ws.p.fill(y_full);
    compute_rates(&ws.p, dose_rate, mat, out);
}

/// Species balance `S·r` [mol/L/s] of the dynamic species.
pub fn species_rates(rates: &ReactionRates, mat: &SimMatrices, acc: &mut [f64]) {
    acc.fill(0.0);
    for (col, r) in mat.k.stoichiometry.iter().zip(&rates.k) {
        for &(i, c) in col {
            acc[i] += c * r;
        }
    }
    for (col, r) in mat.enzymatic.stoichiometry.iter().zip(&rates.enzymatic) {
        for &(i, c) in col {
            acc[i] += c * r;
        }
    }
    for (sp, r) in mat.radiolytic.species.iter().zip(&rates.radiolytic) {
        if let Some(i) = sp {
            acc[*i] += r;
        }
    }
}

/// Right-hand side: `dydt` [µmol/L/s] of the dynamic state `y` [µmol/L].
pub fn derive(y: &[f64], dose_rate: f64, mat: &SimMatrices, ws: &mut Workspace, dydt: &mut [f64]) {
    to_full_molar(y, mat, &mut ws.y_full);
    ws.p.fill(&ws.y_full);
    compute_rates(&ws.p, dose_rate, mat, &mut ws.rates);
    species_rates(&ws.rates, mat, &mut ws.acc);
    for (j, &yj) in y.iter().enumerate().take(mat.n_dyn) {
        if yj < 0.0 {
            let slope = overshoot(j, mat, ws, false);
            ws.acc[j] += slope * ws.y_full[j];
        }
    }
    for (o, a) in dydt.iter_mut().zip(&ws.acc) {
        *o = a * 1e6;
    }
}

fn k_relaxation(e: &KOvershootEntry, p0: &Prepared, grad: GradBuffers<'_>) -> f64 {
    let law = mass_action(&p0.clip, e.k, &e.slots);
    let law_j = mass_action(&p0.clip, e.k_dj, &e.slots_dj);
    let sw = switch(&p0.sw, &e.negative_without_j);
    let effective = if e.flag { law / EPS } else { law_j };
    if let Some((g, d)) = grad {
        d.clear();
        let a = e.weight * sw;
        if e.flag {
            mass_action_derivatives(p0, e.k, &e.slots, d, a / EPS);
        } else {
            mass_action_derivatives(p0, e.k_dj, &e.slots_dj, d, a);
        }
        switch_derivatives(p0, &e.negative_without_j, d, e.weight * effective);
        for &(s, v) in d.iter() {
            g[s] += v;
        }
    }
    e.weight * effective * sw
}

fn ez_relaxation(e: &EzOvershootEntry, p0: &Prepared, grad: GradBuffers<'_>) -> f64 {
    let en = p0.clip[e.enzyme];
    let su = p0.clip[e.substrate];
    let den = e.km + su;
    let d_enzyme = e.k * su / den;
    let d_substrate = e.k * en * e.km / (den * den);
    let law_j = (if e.is_enzyme { d_enzyme } else { 0.0 })
        + (if e.is_substrate { d_substrate } else { 0.0 });
    let law = e.k * en * su / den;
    let sw = switch(&p0.sw, &e.negative_without_j);
    let effective = if e.flag { law / EPS } else { law_j };
    if let Some((g, d)) = grad {
        d.clear();
        let a = e.weight * sw;
        if e.flag {
            d.push((e.enzyme, d_enzyme * a / EPS));
            d.push((e.substrate, d_substrate * a / EPS));
        } else {
            let d_es = e.k * e.km / (den * den);
            let d_ss = -2.0 * e.k * en * e.km / (den * den * den);
            d.push((e.enzyme, (if e.is_substrate { d_es } else { 0.0 }) * a));
            d.push((
                e.substrate,
                ((if e.is_enzyme { d_es } else { 0.0 })
                    + (if e.is_substrate { d_ss } else { 0.0 }))
                    * a,
            ));
        }
        d.retain(|(s, _)| p0.raw[*s] >= 0.0);
        switch_derivatives(p0, &e.negative_without_j, d, e.weight * effective);
        for &(s, v) in d.iter() {
            g[s] += v;
        }
    }
    e.weight * effective * sw
}

/// Slope `J_jj(0⁺)` of species `j` (and its gradient into `ws.grad` when
/// `with_grad`), with `ws.p` prepared for the current state.
fn overshoot(j: usize, mat: &SimMatrices, ws: &mut Workspace, with_grad: bool) -> f64 {
    ws.p0.copy_from(&ws.p);
    ws.p0.set(j, 0.0);
    if with_grad {
        ws.grad.fill(0.0);
    }
    let mut slope = 0.0;
    let Workspace {
        p0, grad, deriv, ..
    } = ws;
    for e in &mat.overshoot.k[j] {
        slope += k_relaxation(e, p0, with_grad.then_some((&mut *grad, &mut *deriv)));
    }
    for e in &mat.overshoot.ez[j] {
        slope += ez_relaxation(e, p0, with_grad.then_some((&mut *grad, &mut *deriv)));
    }
    if with_grad {
        grad[j] = 0.0;
    }
    slope
}

/// Dense analytic Jacobian of [`derive`], column-major
/// `pd[i + j*n_dyn] = ∂f_i/∂y_j` (rsoda layout).
pub fn jacobian(y: &[f64], mat: &SimMatrices, ws: &mut Workspace, pd: &mut [f64]) {
    let n = mat.n_dyn;
    to_full_molar(y, mat, &mut ws.y_full);
    ws.p.fill(&ws.y_full);
    pd.fill(0.0);
    let k = &mat.k;
    for r in 0..k.k_values.len() {
        let col = &k.stoichiometry[r];
        if col.is_empty() {
            continue;
        }
        let neg = &k.negative_index[r];
        let law = mass_action(&ws.p.clip, k.k_values[r], &k.reactant_index[r]);
        let sw = switch(&ws.p.sw, neg);
        ws.deriv.clear();
        mass_action_derivatives(
            &ws.p,
            k.k_values[r],
            &k.reactant_index[r],
            &mut ws.deriv,
            sw,
        );
        switch_derivatives(&ws.p, neg, &mut ws.deriv, law);
        scatter(col, &ws.deriv, n, pd);
    }
    let ez = &mat.enzymatic;
    for r in 0..ez.k_values.len() {
        let col = &ez.stoichiometry[r];
        if col.is_empty() {
            continue;
        }
        let (ei, si) = (ez.enzyme_index[r], ez.substrate_index[r]);
        let (kv, km) = (ez.k_values[r], ez.k_michaelis[r]);
        let (en, su) = (ws.p.clip[ei], ws.p.clip[si]);
        let den = km + su;
        let d_enzyme = if ws.p.raw[ei] < 0.0 {
            0.0
        } else {
            kv * su / den
        };
        let d_substrate = if ws.p.raw[si] < 0.0 {
            0.0
        } else {
            kv * en * km / (den * den)
        };
        let neg = &ez.negative_index[r];
        let law = kv * en * su / den;
        let sw = switch(&ws.p.sw, neg);
        ws.deriv.clear();
        ws.deriv.push((ei, d_enzyme * sw));
        ws.deriv.push((si, d_substrate * sw));
        switch_derivatives(&ws.p, neg, &mut ws.deriv, law);
        scatter(col, &ws.deriv, n, pd);
    }
    for j in 0..n {
        if y[j] < 0.0 {
            let slope = overshoot(j, mat, ws, true);
            let yj = ws.y_full[j];
            for s in 0..n {
                pd[j + s * n] += yj * ws.grad[s];
            }
            pd[j + j * n] += slope;
        }
    }
}

/// Where one reaction writes into the sparse Jacobian: its candidate
/// columns (the dynamic species its rate depends on) and, for each of them,
/// the value positions of its stoichiometry rows,
/// `pos[c * rows + r]` for column `cols[c]` and row `stoichiometry[r]`.
#[derive(Debug, Clone, PartialEq)]
struct Slots {
    cols: Vec<usize>,
    pos: Vec<usize>,
}

impl Slots {
    fn new(rows: &[usize], cols: Vec<usize>, pattern: &SparsityPattern) -> Self {
        let mut pos = Vec::with_capacity(cols.len() * rows.len());
        for &c in &cols {
            for &r in rows {
                pos.push(pattern.find(r, c).expect("entry in the Jacobian pattern"));
            }
        }
        Self { cols, pos }
    }
}

/// Sorted dynamic species of `ids` (below `n`), without duplicates.
fn dyn_columns(ids: impl IntoIterator<Item = usize>, n: usize) -> Vec<usize> {
    let mut v: Vec<usize> = ids.into_iter().filter(|&s| s < n).collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Sparsity pattern of the Jacobian of [`derive`] and the value positions
/// [`jacobian_sparse`] writes to, built once per [`SimMatrices`].
///
/// The pattern is the union of: for each reaction, its stoichiometry rows
/// times the dynamic species its rate depends on (reactants or enzyme and
/// substrate, plus the "negative" products of the depletion switch); for
/// each dynamic species `j`, the columns its overshoot relaxation depends on
/// (row `j`); and the diagonal.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseJacobian {
    pattern: SparsityPattern,
    k: Vec<Slots>,
    ez: Vec<Slots>,
    overshoot: Vec<Slots>,
}

impl SparseJacobian {
    /// # Panics
    ///
    /// If `mat` has no dynamic species.
    pub fn new(mat: &SimMatrices) -> Self {
        let n = mat.n_dyn;
        let (k, ez) = (&mat.k, &mat.enzymatic);
        let rows = |col: &[(usize, f64)]| col.iter().map(|&(i, _)| i).collect::<Vec<_>>();
        let k_cols: Vec<Vec<usize>> = (0..k.k_values.len())
            .map(|r| {
                dyn_columns(
                    k.reactant_index[r]
                        .iter()
                        .chain(&k.negative_index[r])
                        .copied(),
                    n,
                )
            })
            .collect();
        let ez_cols: Vec<Vec<usize>> = (0..ez.k_values.len())
            .map(|r| {
                let es = [ez.enzyme_index[r], ez.substrate_index[r]];
                dyn_columns(
                    es.into_iter().chain(ez.negative_index[r].iter().copied()),
                    n,
                )
            })
            .collect();
        let over_cols: Vec<Vec<usize>> = (0..n)
            .map(|j| {
                let mut ids = Vec::new();
                for e in &mat.overshoot.k[j] {
                    ids.extend(
                        e.slots
                            .iter()
                            .chain(&e.slots_dj)
                            .chain(&e.negative_without_j),
                    );
                }
                for e in &mat.overshoot.ez[j] {
                    ids.extend([e.enzyme, e.substrate]);
                    ids.extend(&e.negative_without_j);
                }
                dyn_columns(ids, n)
            })
            .collect();
        let mut entries = Vec::new();
        for (cols, st) in k_cols
            .iter()
            .zip(&k.stoichiometry)
            .chain(ez_cols.iter().zip(&ez.stoichiometry))
        {
            for &c in cols {
                entries.extend(st.iter().map(|&(i, _)| (i, c)));
            }
        }
        for (j, cols) in over_cols.iter().enumerate() {
            entries.extend(cols.iter().map(|&c| (j, c)));
        }
        let pattern = SparsityPattern::from_entries(n, entries)
            .expect("Jacobian pattern of the dynamic species");
        let k_slots = k_cols
            .into_iter()
            .zip(&k.stoichiometry)
            .map(|(c, st)| Slots::new(&rows(st), c, &pattern))
            .collect();
        let ez_slots = ez_cols
            .into_iter()
            .zip(&ez.stoichiometry)
            .map(|(c, st)| Slots::new(&rows(st), c, &pattern))
            .collect();
        let over_slots = over_cols
            .into_iter()
            .enumerate()
            .map(|(j, c)| Slots::new(&[j], c, &pattern))
            .collect();
        Self {
            pattern,
            k: k_slots,
            ez: ez_slots,
            overshoot: over_slots,
        }
    }

    /// The pattern (CSC, diagonal included) the values of
    /// [`jacobian_sparse`] refer to.
    pub fn pattern(&self) -> &SparsityPattern {
        &self.pattern
    }
}

/// [`jacobian`] on the entries of `sj.pattern()`: `values[k] = ∂f_i/∂y_j`
/// for the `k`-th stored entry `(i, j)` (rsoda's
/// [`rsoda::OdeSystem::jac_sparse`] layout). Same arithmetic as the dense
/// Jacobian, no allocation.
pub fn jacobian_sparse(
    y: &[f64],
    mat: &SimMatrices,
    sj: &SparseJacobian,
    ws: &mut Workspace,
    values: &mut [f64],
) {
    let n = mat.n_dyn;
    to_full_molar(y, mat, &mut ws.y_full);
    ws.p.fill(&ws.y_full);
    values.fill(0.0);
    let k = &mat.k;
    for r in 0..k.k_values.len() {
        let col = &k.stoichiometry[r];
        if col.is_empty() {
            continue;
        }
        let neg = &k.negative_index[r];
        let law = mass_action(&ws.p.clip, k.k_values[r], &k.reactant_index[r]);
        let sw = switch(&ws.p.sw, neg);
        ws.deriv.clear();
        mass_action_derivatives(
            &ws.p,
            k.k_values[r],
            &k.reactant_index[r],
            &mut ws.deriv,
            sw,
        );
        switch_derivatives(&ws.p, neg, &mut ws.deriv, law);
        scatter_sparse(col, &ws.deriv, n, &sj.k[r], values);
    }
    let ez = &mat.enzymatic;
    for r in 0..ez.k_values.len() {
        let col = &ez.stoichiometry[r];
        if col.is_empty() {
            continue;
        }
        let (ei, si) = (ez.enzyme_index[r], ez.substrate_index[r]);
        let (kv, km) = (ez.k_values[r], ez.k_michaelis[r]);
        let (en, su) = (ws.p.clip[ei], ws.p.clip[si]);
        let den = km + su;
        let d_enzyme = if ws.p.raw[ei] < 0.0 {
            0.0
        } else {
            kv * su / den
        };
        let d_substrate = if ws.p.raw[si] < 0.0 {
            0.0
        } else {
            kv * en * km / (den * den)
        };
        let neg = &ez.negative_index[r];
        let law = kv * en * su / den;
        let sw = switch(&ws.p.sw, neg);
        ws.deriv.clear();
        ws.deriv.push((ei, d_enzyme * sw));
        ws.deriv.push((si, d_substrate * sw));
        switch_derivatives(&ws.p, neg, &mut ws.deriv, law);
        scatter_sparse(col, &ws.deriv, n, &sj.ez[r], values);
    }
    let diag = sj.pattern.diag();
    for j in 0..n {
        if y[j] < 0.0 {
            let slope = overshoot(j, mat, ws, true);
            let yj = ws.y_full[j];
            let row = &sj.overshoot[j];
            for (&s, &k) in row.cols.iter().zip(&row.pos) {
                values[k] += yj * ws.grad[s];
            }
            values[diag[j]] += slope;
        }
    }
}

/// [`scatter`] into the value positions of one reaction.
fn scatter_sparse(
    col: &[(usize, f64)],
    deriv: &[(usize, f64)],
    n: usize,
    slots: &Slots,
    values: &mut [f64],
) {
    let rows = col.len();
    for &(s, d) in deriv {
        if s >= n {
            continue;
        }
        let c = slots
            .cols
            .iter()
            .position(|&x| x == s)
            .expect("derivative column in the reaction's slots");
        let pos = &slots.pos[c * rows..(c + 1) * rows];
        for (&(_, coef), &k) in col.iter().zip(pos) {
            values[k] += coef * d;
        }
    }
}

fn scatter(col: &[(usize, f64)], deriv: &[(usize, f64)], n: usize, pd: &mut [f64]) {
    for &(s, d) in deriv {
        if s >= n {
            continue;
        }
        for &(i, c) in col {
            pd[i + s * n] += c * d;
        }
    }
}
