//! Reaction-diffusion master equation on the voxel grid (ADR 0003).
//!
//! Every voxel holds integer molecule counts. Two kinds of events compete:
//! a reaction inside a voxel, with the propensities of [`SsaModel`] for the
//! voxel volume, and a hop of one molecule to a face neighbour, with
//! propensity `D / h^2` per molecule and face (zero-flux boundaries: a
//! boundary voxel has fewer faces). The mean of the hop process obeys the
//! same discrete diffusion operator as the deterministic voxel path.
//!
//! The simulation is exact (next-subvolume method): each voxel has the time
//! of its next event in a priority queue; the earliest voxel fires one event,
//! chosen in proportion to the rates, and only the voxels whose counts
//! changed draw a new time.

use crate::error::{RadkinError, Result};
use crate::matrices::SimMatrices;
use crate::ssa::model::SsaModel;
use crate::ssa::rng::SsaRng;
use crate::voxel::coarsen::coarse_grid;
use crate::voxel::grid::VoxelGrid;

/// Next event time of every voxel: an indexed binary min-heap.
#[derive(Debug, Clone, Default)]
pub(crate) struct EventQueue {
    key: Vec<f64>,
    heap: Vec<u32>,
    pos: Vec<u32>,
}

impl EventQueue {
    pub(crate) fn build(keys: Vec<f64>) -> Self {
        let n = keys.len();
        let mut q = Self {
            key: keys,
            heap: (0..n as u32).collect(),
            pos: (0..n as u32).collect(),
        };
        for i in (0..n / 2).rev() {
            q.sift_down(i);
        }
        q
    }

    /// The voxel with the earliest event and its time.
    pub(crate) fn top(&self) -> (usize, f64) {
        let v = self.heap[0] as usize;
        (v, self.key[v])
    }

    pub(crate) fn set(&mut self, v: usize, key: f64) {
        let old = self.key[v];
        self.key[v] = key;
        let i = self.pos[v] as usize;
        if key < old {
            self.sift_up(i);
        } else {
            self.sift_down(i);
        }
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.heap.swap(i, j);
        self.pos[self.heap[i] as usize] = i as u32;
        self.pos[self.heap[j] as usize] = j as u32;
    }

    fn sift_up(&mut self, mut i: usize) {
        while i > 0 {
            let parent = (i - 1) / 2;
            if self.key[self.heap[i] as usize] < self.key[self.heap[parent] as usize] {
                self.swap(i, parent);
                i = parent;
            } else {
                break;
            }
        }
    }

    fn sift_down(&mut self, mut i: usize) {
        let n = self.heap.len();
        loop {
            let (l, r) = (2 * i + 1, 2 * i + 2);
            let mut m = i;
            if l < n && self.key[self.heap[l] as usize] < self.key[self.heap[m] as usize] {
                m = l;
            }
            if r < n && self.key[self.heap[r] as usize] < self.key[self.heap[m] as usize] {
                m = r;
            }
            if m == i {
                break;
            }
            self.swap(i, m);
            i = m;
        }
    }
}

/// The grid and the propensity model of one coarsening level.
#[derive(Debug, Clone)]
pub(crate) struct Level {
    pub(crate) grid: VoxelGrid,
    pub(crate) model: SsaModel,
    /// `1 / h^2` per axis [1/m^2].
    pub(crate) inv_h2: [f64; 3],
}

/// `C_3 / 6` with `C_3 = 1.5164` (Hellander, Hellander and Petzold 2012,
/// Theorem 1): the critical spacing of a pair is `h* = C_3 kappa / (6 D)`.
const CRITICAL_FACTOR: f64 = 1.5164 / 6.0;

/// Critical spacing `h*` [m] of every mass-action reaction with two dynamic
/// reactant slots and a diffusing pair, as `(reaction, h*)`.
///
/// `kappa` [m^3/s] is the pair rate constant the RDME sees: `c N_A V`
/// (`k` times any fixed concentration) per molecule, doubled for `2 A`,
/// whose propensity `c n (n - 1)` is `2 c` per pair. `D` is the summed
/// diffusion coefficient of the pair (`2 D_A` for `2 A`).
fn critical_spacings(model: &SsaModel, diffusion: &[f64]) -> Vec<(usize, f64)> {
    use crate::consts::AVOGADRO;
    use crate::ssa::model::Law;
    let mut out = Vec::new();
    for (j, r) in model.reactions.iter().enumerate() {
        let Law::MassAction { c, reactants } = &r.law else {
            continue;
        };
        if r.order != 2 {
            continue;
        }
        let d: f64 = reactants
            .iter()
            .map(|(s, mult)| diffusion[*s] * *mult as f64)
            .sum();
        if d <= 0.0 {
            continue;
        }
        let pairs = if reactants.len() == 1 { 2.0 } else { 1.0 };
        let kappa = pairs * c * model.nav / (1e3 * AVOGADRO);
        out.push((j, CRITICAL_FACTOR * kappa / d));
    }
    out
}

/// Whether the three spacings of `grid` are equal (relative 1e-12).
pub fn is_cubic(grid: &VoxelGrid) -> bool {
    let h = grid.spacing[0];
    grid.spacing.iter().all(|s| (s - h).abs() <= 1e-12 * h)
}

/// Everything that does not change during a run, shared by the trajectories.
#[derive(Debug, Clone)]
pub struct RdmeSetup {
    pub(crate) levels: Vec<Level>,
    /// Diffusion coefficient [m^2/s] of every dynamic species.
    pub(crate) diffusion: Vec<f64>,
    pub(crate) n_sp: usize,
    /// `(reaction, h*)` of every bimolecular reaction with a diffusing pair.
    critical: Vec<(usize, f64)>,
    /// Whether the propensities carry the mesh correction.
    corrected: bool,
}

impl RdmeSetup {
    /// Setup for `mat` on `grid`, with `max_levels` coarsening levels above
    /// the fine grid (every dimension must be divisible by `2^max_levels`).
    ///
    /// With `correct`, every level of spacing `h` scales the propensity of
    /// each reaction of critical spacing `h* < h` by `1 / (1 - h*/h)`
    /// (mesh correction, Hellander and Hellander 2020, eq. 6 with
    /// `1/k = 1/k_r + 1/(4 pi sigma D)`); the grid must be cubic. Reactions
    /// with `h* >= h` keep the macroscopic rate on that level.
    pub fn new(
        mat: &SimMatrices,
        grid: &VoxelGrid,
        diffusion: Vec<f64>,
        max_levels: u32,
        correct: bool,
    ) -> Result<Self> {
        if diffusion.len() != mat.n_dyn {
            return Err(RadkinError::Value(format!(
                "rdme: {} diffusion coefficients for {} dynamic species",
                diffusion.len(),
                mat.n_dyn
            )));
        }
        if diffusion.iter().any(|d| !(d.is_finite() && *d >= 0.0)) {
            return Err(RadkinError::Value(
                "rdme: diffusion coefficients must be finite and >= 0".into(),
            ));
        }
        if correct && !is_cubic(grid) {
            return Err(RadkinError::Value(format!(
                "rdme: the mesh correction needs a cubic grid, got spacing {:?} m",
                grid.spacing
            )));
        }
        // `h*` does not depend on the volume: take it from the fine model.
        let critical = critical_spacings(&SsaModel::new(mat, grid.voxel_volume_l())?, &diffusion);
        let mut levels = Vec::with_capacity(max_levels as usize + 1);
        let mut g = grid.clone();
        for l in 0..=max_levels {
            if l > 0 {
                g = coarse_grid(&g)?;
            }
            let mut model = SsaModel::new(mat, g.voxel_volume_l())?;
            if correct {
                mesh_correct(&mut model, &critical, g.spacing[0]);
            }
            levels.push(Level {
                model,
                inv_h2: g.spacing.map(|h| 1.0 / (h * h)),
                grid: g.clone(),
            });
        }
        Ok(Self {
            levels,
            diffusion,
            n_sp: mat.n_dyn,
            critical,
            corrected: correct,
        })
    }

    /// Whether the propensities carry the mesh correction.
    pub fn corrected(&self) -> bool {
        self.corrected
    }

    /// `(reaction, h*)` of every bimolecular reaction with a diffusing pair:
    /// its critical spacing [m].
    pub fn critical_spacings(&self) -> &[(usize, f64)] {
        &self.critical
    }

    /// The fine grid.
    pub fn grid(&self) -> &VoxelGrid {
        &self.levels[0].grid
    }

    /// The propensity model of a fine voxel.
    pub fn model(&self) -> &SsaModel {
        &self.levels[0].model
    }

    /// Coarsening levels above the fine grid.
    pub fn max_levels(&self) -> u32 {
        self.levels.len() as u32 - 1
    }

    /// Smallest voxel spacing below which a bimolecular reaction can no
    /// longer be represented by the RDME, and the index of that reaction:
    /// the largest critical spacing `h* = 0.2527 kappa / D`, `kappa` the pair
    /// rate constant and `D` the sum of the two diffusion coefficients
    /// (Hellander, Hellander and Petzold, Phys. Rev. E 85, 042901, 2012).
    /// `None` when no reaction has two diffusing reactants.
    pub fn validity_bound(&self) -> Option<(f64, usize)> {
        let mut worst: Option<(f64, usize)> = None;
        for &(j, h) in &self.critical {
            if worst.is_none_or(|(w, _)| h > w) {
                worst = Some((h, j));
            }
        }
        worst
    }
}

/// Scale the propensity of every reaction of `critical` with `h* < h` by
/// `1 / (1 - h*/h)`.
fn mesh_correct(model: &mut SsaModel, critical: &[(usize, f64)], h: f64) {
    use crate::ssa::model::Law;
    for &(j, h_star) in critical {
        if h_star >= h {
            continue;
        }
        if let Law::MassAction { c, .. } = &mut model.reactions[j].law {
            *c /= 1.0 - h_star / h;
        }
    }
}

/// One trajectory of the RDME.
#[derive(Debug, Clone)]
pub struct Rdme<'a> {
    pub(crate) setup: &'a RdmeSetup,
    /// Current coarsening level (0: the fine grid).
    pub(crate) level: usize,
    pub t: f64,
    /// Counts of the current level, voxel-major: `counts[v * n_sp + s]`.
    pub(crate) counts: Vec<i64>,
    /// Total reaction propensity of every voxel.
    react: Vec<f64>,
    /// `sum_s D_s n_s` of every voxel.
    dsum: Vec<f64>,
    queue: EventQueue,
    pub(crate) rng: SsaRng,
    /// Scratch: propensities of one voxel.
    a: Vec<f64>,
    pub n_reactions: u64,
    pub n_hops: u64,
    /// Events selected in a voxel whose rate had become zero (round-off of
    /// the incremental rate updates); they change nothing.
    pub n_null: u64,
    /// Threshold of the dispersion test; `None`: coarsening off.
    pub(crate) z_max: Option<f64>,
    /// `(t, level)` after every change of coarsening level.
    pub level_history: Vec<(f64, u32)>,
    /// Events fired since the last coarsening test.
    events_since_test: u64,
    /// Tau-leaping state; `None`: exact.
    pub(crate) tau: Option<Box<crate::ssa::rdme_tau::RdmeTau>>,
    /// The event queue does not match the counts (set by a leap).
    pub(crate) queue_stale: bool,
}

/// Stride of every axis in the flat voxel index.
pub(crate) fn strides(grid: &VoxelGrid) -> [usize; 3] {
    [1, grid.dims[0], grid.dims[0] * grid.dims[1]]
}

/// `sum over the faces of voxel v of 1 / h^2`.
pub(crate) fn face_weight(grid: &VoxelGrid, inv_h2: &[f64; 3], v: usize) -> f64 {
    let c = grid.coords(v);
    let mut w = 0.0;
    for a in 0..3 {
        if c[a] > 0 {
            w += inv_h2[a];
        }
        if c[a] + 1 < grid.dims[a] {
            w += inv_h2[a];
        }
    }
    w
}

impl<'a> Rdme<'a> {
    /// An empty grid at time `t`.
    pub fn new(setup: &'a RdmeSetup, t: f64, rng: SsaRng) -> Self {
        let n_vox = setup.grid().n_voxels();
        let mut r = Self {
            setup,
            level: 0,
            t,
            counts: vec![0; n_vox * setup.n_sp],
            react: Vec::new(),
            dsum: Vec::new(),
            queue: EventQueue::default(),
            rng,
            a: vec![0.0; setup.model().reactions.len()],
            n_reactions: 0,
            n_hops: 0,
            n_null: 0,
            z_max: None,
            level_history: Vec::new(),
            events_since_test: 0,
            tau: None,
            queue_stale: false,
        };
        r.rebuild();
        r
    }

    /// Grid of the current level.
    pub fn grid(&self) -> &VoxelGrid {
        &self.setup.levels[self.level].grid
    }

    /// Current coarsening level.
    pub fn level(&self) -> u32 {
        self.level as u32
    }

    /// Reaction and hop events fired.
    pub fn n_events(&self) -> u64 {
        self.n_reactions + self.n_hops
    }

    /// Recompute every voxel rate and draw every next event time.
    pub(crate) fn rebuild(&mut self) {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let (n_sp, n_vox) = (self.setup.n_sp, lv.grid.n_voxels());
        debug_assert_eq!(self.counts.len(), n_sp * n_vox);
        self.react.clear();
        self.dsum.clear();
        let mut keys = Vec::with_capacity(n_vox);
        for v in 0..n_vox {
            let n = &self.counts[v * n_sp..(v + 1) * n_sp];
            let react: f64 = (0..lv.model.reactions.len())
                .map(|j| lv.model.propensity(j, n))
                .sum();
            let dsum: f64 = n
                .iter()
                .zip(&self.setup.diffusion)
                .map(|(c, d)| *c as f64 * d)
                .sum();
            let rate = react + dsum * face_weight(&lv.grid, &lv.inv_h2, v);
            keys.push(if rate > 0.0 {
                self.t + self.rng.exp1() / rate
            } else {
                f64::INFINITY
            });
            self.react.push(react);
            self.dsum.push(dsum);
        }
        self.queue = EventQueue::build(keys);
        self.queue_stale = false;
    }

    /// Set every count to a Poisson draw of mean `c N_A V` for the uniform
    /// concentrations `y` [µmol/L] (a uniform solute is Poisson distributed
    /// over the voxels).
    pub fn fill_poisson(&mut self, y: &[f64]) {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let k = lv.model.molecules_per_umol();
        let n_sp = self.setup.n_sp;
        assert_eq!(y.len(), n_sp, "rdme: state length");
        for v in 0..lv.grid.n_voxels() {
            for (s, c) in y.iter().enumerate() {
                self.counts[v * n_sp + s] = self.rng.poisson(c * k);
            }
        }
        self.rebuild();
    }

    /// Add `(fine voxel, species, molecules)` entries. A coarsened grid is
    /// refined to the fine grid first.
    pub fn inject(&mut self, molecules: &[(usize, usize, i64)]) {
        self.refine_to_fine();
        self.events_since_test = 0;
        let n_sp = self.setup.n_sp;
        for &(v, s, n) in molecules {
            self.counts[v * n_sp + s] += n;
        }
        self.rebuild();
    }

    /// Molecules of every species in the domain.
    pub fn totals(&self) -> Vec<i64> {
        let n_sp = self.setup.n_sp;
        let mut out = vec![0; n_sp];
        for voxel in self.counts.chunks_exact(n_sp) {
            for (o, c) in out.iter_mut().zip(voxel) {
                *o += c;
            }
        }
        out
    }

    /// Counts of voxel `v` of the current level.
    pub fn voxel(&self, v: usize) -> &[i64] {
        let n_sp = self.setup.n_sp;
        &self.counts[v * n_sp..(v + 1) * n_sp]
    }

    /// The field of every species on the fine grid, species-major
    /// (`out[s * n_fine + v]`), in units of `1 / 8^max_levels` molecule: a
    /// voxel of level `L` spreads its count evenly over its `8^L` fine
    /// voxels, which stays an integer in these units.
    pub fn fine_field(&self) -> Vec<i64> {
        let n_sp = self.setup.n_sp;
        let fine = self.setup.grid();
        let n_fine = fine.n_voxels();
        let max = self.setup.max_levels() as usize;
        let scale = 8i64.pow((max - self.level) as u32);
        let shift = self.level;
        let cur = self.grid();
        let mut out = vec![0; n_sp * n_fine];
        for v in 0..n_fine {
            let c = fine.coords(v);
            let cv = cur.index(c[0] >> shift, c[1] >> shift, c[2] >> shift);
            for s in 0..n_sp {
                out[s * n_fine + v] = self.counts[cv * n_sp + s] * scale;
            }
        }
        out
    }

    /// New next-event time of voxel `v` from its current rates.
    fn reschedule(&mut self, v: usize) {
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let rate = self.react[v] + self.dsum[v] * face_weight(&lv.grid, &lv.inv_h2, v);
        let key = if rate > 0.0 {
            self.t + self.rng.exp1() / rate
        } else {
            f64::INFINITY
        };
        self.queue.set(v, key);
    }

    /// Exact `sum_s D_s n_s` of voxel `v`.
    fn refresh_dsum(&mut self, v: usize) {
        let n_sp = self.setup.n_sp;
        self.dsum[v] = self.counts[v * n_sp..(v + 1) * n_sp]
            .iter()
            .zip(&self.setup.diffusion)
            .map(|(c, d)| *c as f64 * d)
            .sum();
    }

    /// Sum of the propensities of the reactions that read species `s` in voxel `v`.
    fn readers_sum(&self, v: usize, s: usize) -> f64 {
        let setup = self.setup;
        let model = &setup.levels[self.level].model;
        let n_sp = self.setup.n_sp;
        let n = &self.counts[v * n_sp..(v + 1) * n_sp];
        model.readers[s]
            .iter()
            .map(|j| model.propensity(*j, n))
            .sum()
    }

    /// Fire a reaction in voxel `v`; `frac` in `[0, 1)` selects it.
    fn fire_reaction(&mut self, v: usize, frac: f64) {
        let setup = self.setup;
        let model = &setup.levels[self.level].model;
        let n_sp = self.setup.n_sp;
        let n = &self.counts[v * n_sp..(v + 1) * n_sp];
        let mut total = 0.0;
        for (j, a) in self.a.iter_mut().enumerate() {
            *a = model.propensity(j, n);
            total += *a;
        }
        if total <= 0.0 {
            // The incremental rate had drifted above an exact zero.
            self.react[v] = 0.0;
            self.n_null += 1;
            return;
        }
        let target = frac * total;
        let (mut acc, mut chosen) = (0.0, 0);
        for (j, a) in self.a.iter().enumerate() {
            if *a > 0.0 {
                acc += a;
                chosen = j;
                if target < acc {
                    break;
                }
            }
        }
        for &(s, d) in &model.reactions[chosen].delta {
            self.counts[v * n_sp + s] += d;
            debug_assert!(self.counts[v * n_sp + s] >= 0);
        }
        let n = &self.counts[v * n_sp..(v + 1) * n_sp];
        self.react[v] = (0..model.reactions.len())
            .map(|j| model.propensity(j, n))
            .sum();
        self.refresh_dsum(v);
        self.n_reactions += 1;
    }

    /// Fire a hop out of voxel `v`; `x` in `[0, dsum)` selects the species.
    /// Returns the destination voxel, or `None` for a null event.
    fn fire_hop(&mut self, v: usize, x: f64) -> Option<usize> {
        let n_sp = self.setup.n_sp;
        let (mut acc, mut species) = (0.0, None);
        for s in 0..n_sp {
            let n = self.counts[v * n_sp + s];
            if n > 0 && self.setup.diffusion[s] > 0.0 {
                acc += n as f64 * self.setup.diffusion[s];
                species = Some(s);
                if x < acc {
                    break;
                }
            }
        }
        let Some(s) = species else {
            self.dsum[v] = 0.0;
            self.n_null += 1;
            return None;
        };
        // Destination: a face of the voxel, in proportion to 1 / h^2.
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let c = lv.grid.coords(v);
        let st = strides(&lv.grid);
        let target = self.rng.uniform() * face_weight(&lv.grid, &lv.inv_h2, v);
        let (mut acc, mut dest) = (0.0, v);
        'faces: for a in 0..3 {
            if c[a] > 0 {
                acc += lv.inv_h2[a];
                dest = v - st[a];
                if target < acc {
                    break 'faces;
                }
            }
            if c[a] + 1 < lv.grid.dims[a] {
                acc += lv.inv_h2[a];
                dest = v + st[a];
                if target < acc {
                    break 'faces;
                }
            }
        }
        if dest == v {
            // A single-voxel grid has no face.
            self.n_null += 1;
            return None;
        }
        let before = (self.readers_sum(v, s), self.readers_sum(dest, s));
        self.counts[v * n_sp + s] -= 1;
        self.counts[dest * n_sp + s] += 1;
        let after = (self.readers_sum(v, s), self.readers_sum(dest, s));
        self.react[v] = (self.react[v] + (after.0 - before.0)).max(0.0);
        self.react[dest] = (self.react[dest] + (after.1 - before.1)).max(0.0);
        self.refresh_dsum(v);
        self.refresh_dsum(dest);
        self.n_hops += 1;
        Some(dest)
    }

    /// Fire the next event if it happens before `t_stop`. Returns whether
    /// one did.
    pub(crate) fn step(&mut self, t_stop: f64) -> bool {
        let (v, t_next) = self.queue.top();
        if t_next >= t_stop {
            return false;
        }
        self.t = t_next;
        let setup = self.setup;
        let lv = &setup.levels[self.level];
        let hop = self.dsum[v] * face_weight(&lv.grid, &lv.inv_h2, v);
        let u = self.rng.uniform() * (self.react[v] + hop);
        if u < hop {
            let x = u / hop * self.dsum[v];
            if let Some(dest) = self.fire_hop(v, x) {
                self.reschedule(dest);
            }
        } else {
            let frac = ((u - hop) / self.react[v]).clamp(0.0, 1.0 - f64::EPSILON);
            self.fire_reaction(v, frac);
        }
        self.reschedule(v);
        true
    }

    /// Events between two coarsening tests: the test reads every count once,
    /// which costs about as much as a few hundred events.
    fn test_interval(&self) -> u64 {
        8 * self.grid().n_voxels() as u64
    }

    /// Run the coarsening test if enough events fired since the last one.
    pub(crate) fn count_events(&mut self, fired: u64) {
        if self.z_max.is_none() {
            return;
        }
        self.events_since_test += fired;
        if self.events_since_test >= self.test_interval() {
            self.events_since_test = 0;
            self.coarsen_while_uniform();
        }
    }

    /// Advance to `t_stop`: exact events, or tau-leaping where it pays
    /// when it is switched on.
    pub fn advance_to(&mut self, t_stop: f64) {
        if self.tau.is_some() {
            self.advance_tau(t_stop);
            return;
        }
        while self.step(t_stop) {
            self.count_events(1);
        }
        if t_stop > self.t {
            self.t = t_stop;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;

    #[test]
    fn event_queue_orders_and_updates() {
        let mut q = EventQueue::build(vec![5.0, 1.0, 3.0, f64::INFINITY, 2.0]);
        assert_eq!(q.top(), (1, 1.0));
        q.set(1, 10.0);
        assert_eq!(q.top(), (4, 2.0));
        q.set(3, 0.5);
        assert_eq!(q.top(), (3, 0.5));
        q.set(3, f64::INFINITY);
        q.set(4, f64::INFINITY);
        assert_eq!(q.top(), (2, 3.0));
        // Drain in order.
        let mut seen = Vec::new();
        loop {
            let (v, t) = q.top();
            if t.is_infinite() {
                break;
            }
            seen.push(t);
            q.set(v, f64::INFINITY);
        }
        assert_eq!(seen, [3.0, 5.0, 10.0]);
    }

    fn setup(toml: &str, dims: [usize; 3], h: f64, d: f64) -> (SimEnv, RdmeSetup) {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new(dims, [h; 3], [0.0; 3]).unwrap();
        let s = RdmeSetup::new(&mat, &grid, vec![d; mat.n_dyn], 0, false).unwrap();
        (env, s)
    }

    /// Setup with `levels` coarsening levels, with or without the mesh
    /// correction.
    fn setup_levels(toml: &str, h: f64, d: f64, levels: u32, correct: bool) -> RdmeSetup {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new([4, 4, 4], [h; 3], [0.0; 3]).unwrap();
        RdmeSetup::new(&mat, &grid, vec![d; mat.n_dyn], levels, correct).unwrap()
    }

    fn c_of(s: &RdmeSetup, level: usize, j: usize) -> f64 {
        match &s.levels[level].model.reactions[j].law {
            crate::ssa::model::Law::MassAction { c, .. } => *c,
            other => panic!("{other:?}"),
        }
    }

    const PAIR: &str =
        "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = 3.0e10";

    /// `h*` of `PAIR` with D = 5e-9 per species.
    fn pair_h_star() -> f64 {
        CRITICAL_FACTOR * 3.0e10 / (1e3 * crate::consts::AVOGADRO) / 1e-8
    }

    const DECAY: &str = r#"
        [[reactions.k_reaction]]
        reaction = "H_r -> H2"
        k_value = 2.0e6
    "#;

    #[test]
    fn hops_conserve_molecules_and_reach_every_voxel() {
        let (env, s) = setup("[reactions.radiolytic]\ne_aq = 2.6", [4, 3, 2], 10e-9, 5e-9);
        let e = env.species.index_of("e_aq").unwrap();
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(1, 0));
        r.inject(&[(0, e, 2400)]);
        // h^2 / D = 20 ns: well mixed after a microsecond.
        r.advance_to(1e-6);
        assert_eq!(r.t, 1e-6);
        assert_eq!(r.totals()[e], 2400);
        assert_eq!(r.n_reactions, 0);
        assert!(r.n_hops > 10_000);
        for v in 0..24 {
            // Uniform: Binomial(2400, 1/24), mean 100, sd 9.8.
            let n = r.voxel(v)[e];
            assert!((n - 100).abs() < 50, "voxel {v}: {n}");
        }
    }

    #[test]
    fn reaction_rate_does_not_depend_on_the_voxels() {
        // First-order decay: the survival probability exp(-k t) is the same
        // whatever the hops do.
        let (env, s) = setup(DECAY, [3, 3, 3], 10e-9, 5e-9);
        let (h, h2) = (
            env.species.index_of("H_r").unwrap(),
            env.species.index_of("H2").unwrap(),
        );
        let (n0, runs) = (500, 200);
        let mut sum = 0.0;
        for stream in 0..runs {
            let mut r = Rdme::new(&s, 0.0, SsaRng::new(3, stream));
            r.inject(&[(13, h, n0)]);
            r.advance_to(3e-7);
            let t = r.totals();
            assert_eq!(t[h] + t[h2], n0);
            sum += t[h] as f64;
        }
        let p = (-2.0e6f64 * 3e-7).exp();
        let mean = sum / runs as f64;
        let sem = (n0 as f64 * p * (1.0 - p) / runs as f64).sqrt();
        assert!((mean - n0 as f64 * p).abs() <= 4.0 * sem, "{mean}");
    }

    #[test]
    fn poisson_background_has_the_right_mean() {
        let (env, s) = setup(DECAY, [6, 6, 6], 100e-9, 0.0);
        let h = env.species.index_of("H_r").unwrap();
        let mut y = vec![0.0; s.n_sp];
        // 5 molecules per voxel on average.
        y[h] = 5.0 / s.model().molecules_per_umol();
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(2, 0));
        r.fill_poisson(&y);
        let total = r.totals()[h] as f64;
        // Poisson(1080): sd 32.9.
        assert!((total - 1080.0).abs() < 4.0 * 32.9, "{total}");
        // No diffusion coefficient: no hop.
        r.advance_to(1e-12);
        assert_eq!(r.n_hops, 0);
    }

    #[test]
    fn fine_field_is_species_major() {
        let (env, s) = setup(DECAY, [2, 1, 1], 10e-9, 0.0);
        let h = env.species.index_of("H_r").unwrap();
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(1, 0));
        r.inject(&[(1, h, 7)]);
        let f = r.fine_field();
        assert_eq!(f[h * 2], 0);
        assert_eq!(f[h * 2 + 1], 7);
        assert_eq!(f.iter().sum::<i64>(), 7);
    }

    #[test]
    fn validity_bound_follows_k_over_d() {
        let (_, s) = setup(PAIR, [2, 2, 2], 10e-9, 5e-9);
        let (h, j) = s.validity_bound().unwrap();
        assert_eq!(j, 0);
        // 0.2527 * 3e10 / (1e3 N_A) / 1e-8 m^2/s = 1.26 nm.
        let expected = pair_h_star();
        assert!((h - expected).abs() < 1e-12 * expected, "{h}");
        assert!((h - 1.259e-9).abs() < 1e-12, "{h}");
        let (_, none) = setup(DECAY, [2, 2, 2], 10e-9, 5e-9);
        assert!(none.validity_bound().is_none());
    }

    #[test]
    fn two_a_uses_the_pair_constant() {
        // radiopyo stores `k = k_value / 2` for `2 A` (rate `k [A]^2`,
        // stoichiometry -2): the propensity c n (n - 1) is 2c per pair, so the
        // pair constant is 2k = k_value, and D = 2 D_A. Same k_value and
        // summed D as PAIR: the same h* (the per-molecule k alone gives half).
        let (_, s) = setup(
            "[[reactions.k_reaction]]\nreaction = \"OH_r + OH_r -> H2O2\"\nk_value = 3.0e10",
            [2, 2, 2],
            10e-9,
            5e-9,
        );
        let c = c_of(&s, 0, 0);
        let k = c * s.model().nav;
        assert!((k - 1.5e10).abs() < 1e-6 * 1.5e10, "{k}");
        let (h, _) = s.validity_bound().unwrap();
        let expected = pair_h_star();
        assert!((h - expected).abs() < 1e-12 * expected, "{h}");
    }

    #[test]
    fn mesh_correction_scales_each_level_with_its_own_spacing() {
        let h = 2e-9;
        let plain = setup_levels(PAIR, h, 5e-9, 1, false);
        let fixed = setup_levels(PAIR, h, 5e-9, 1, true);
        assert!(!plain.corrected() && fixed.corrected());
        let hs = pair_h_star();
        for (level, spacing) in [(0, h), (1, 2.0 * h)] {
            let want = c_of(&plain, level, 0) / (1.0 - hs / spacing);
            let got = c_of(&fixed, level, 0);
            assert!(
                (got - want).abs() < 1e-12 * want,
                "level {level}: {got} vs {want}"
            );
        }
        // The factor is larger on the fine level, closer to the pole.
        let f0 = c_of(&fixed, 0, 0) / c_of(&plain, 0, 0);
        let f1 = c_of(&fixed, 1, 0) / c_of(&plain, 1, 0);
        assert!(f0 > f1 && f1 > 1.0, "{f0} {f1}");
        assert_eq!(plain.critical_spacings(), fixed.critical_spacings());
    }

    #[test]
    fn mesh_correction_leaves_reactions_at_or_below_the_pole() {
        // h = 1 nm < h* = 1.26 nm on the fine level, 2 nm > h* one level up.
        let plain = setup_levels(PAIR, 1e-9, 5e-9, 1, false);
        let fixed = setup_levels(PAIR, 1e-9, 5e-9, 1, true);
        assert_eq!(c_of(&fixed, 0, 0), c_of(&plain, 0, 0));
        assert!(c_of(&fixed, 1, 0) > c_of(&plain, 1, 0));
    }

    #[test]
    fn mesh_correction_skips_first_order_and_still_pairs() {
        let toml = format!("{DECAY}\n{PAIR}");
        // D = 0: no diffusing pair, no correction.
        let still = setup_levels(&toml, 2e-9, 0.0, 0, true);
        assert!(still.critical_spacings().is_empty());
        let plain = setup_levels(&toml, 2e-9, 5e-9, 0, false);
        let fixed = setup_levels(&toml, 2e-9, 5e-9, 0, true);
        // Reaction 0 is the first-order decay.
        assert_eq!(c_of(&fixed, 0, 0), c_of(&plain, 0, 0));
        assert!(c_of(&fixed, 0, 1) > c_of(&plain, 0, 1));
    }

    #[test]
    fn mesh_correction_uses_the_effective_constant_of_a_fixed_reactant() {
        // H2O is fixed (55 M): kappa = k [H2O].
        let toml = "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r + H2O -> OH_minus\"\
                    \nk_value = 1.0e7";
        let s = setup_levels(toml, 10e-9, 5e-9, 0, false);
        let model = &s.levels[0].model;
        assert_eq!(model.reactions[0].order, 2);
        let k_eff = c_of(&s, 0, 0) * model.nav; // k [H2O], L/mol/s
        assert!(k_eff > 1.0e7 * 50.0, "{k_eff}");
        let (_, h_star) = s.critical_spacings()[0];
        let want = CRITICAL_FACTOR * k_eff / (1e3 * crate::consts::AVOGADRO) / 1e-8;
        assert!((h_star - want).abs() < 1e-12 * want, "{h_star}");
    }

    #[test]
    fn mesh_correction_refuses_an_anisotropic_grid() {
        let env = SimEnv::from_config(&Config::from_toml_str(PAIR, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new([2, 2, 2], [2e-9, 2e-9, 4e-9], [0.0; 3]).unwrap();
        assert!(!is_cubic(&grid));
        assert!(RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], 0, true).is_err());
        assert!(RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], 0, false).is_ok());
    }

    /// Mean time until one immobile A (centre voxel) and one B (uniform
    /// start, diffusion `D`) react, over `runs` trajectories, divided by the
    /// macroscopic `V / kappa`; and its standard error. The grid has
    /// `h = 2 h*`.
    fn association_time(correct: bool, runs: u64) -> (f64, f64) {
        let (h, d_b, n) = (6e-9, 5e-9, 11);
        let kappa = (h / 2.0) * d_b / CRITICAL_FACTOR;
        let k = kappa * 1e3 * crate::consts::AVOGADRO;
        let toml = format!(
            "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = {k:e}"
        );
        let env = SimEnv::from_config(&Config::from_toml_str(&toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new([n; 3], [h; 3], [0.0; 3]).unwrap();
        let (a, b) = (
            env.species.index_of("e_aq").unwrap(),
            env.species.index_of("OH_r").unwrap(),
        );
        let mut diffusion = vec![0.0; mat.n_dyn];
        diffusion[b] = d_b;
        let s = RdmeSetup::new(&mat, &grid, diffusion, 0, correct).unwrap();
        let (h_star, _) = s.validity_bound().unwrap();
        assert!((h_star - h / 2.0).abs() < 1e-12 * h);
        let nv = grid.n_voxels();
        let tau = grid.voxel_volume_l() * 1e-3 * nv as f64 / kappa;
        let mut start = SsaRng::new(99, 0);
        let (mut sum, mut sum_sq) = (0.0, 0.0);
        for i in 0..runs {
            let vb = ((start.uniform() * nv as f64) as usize).min(nv - 1);
            let mut r = Rdme::new(&s, 0.0, SsaRng::new(7, i));
            r.inject(&[(nv / 2, a, 1), (vb, b, 1)]);
            while r.n_reactions == 0 {
                assert!(r.step(f64::INFINITY));
            }
            sum += r.t / tau;
            sum_sq += (r.t / tau).powi(2);
        }
        let mean = sum / runs as f64;
        let sem = ((sum_sq / runs as f64 - mean * mean) / runs as f64).sqrt();
        (mean, sem)
    }

    #[test]
    fn mesh_correction_matches_the_association_time() {
        // Hellander, Hellander and Petzold 2012, eq. 1 with N_steps =
        // 1.5164 N: the uncorrected RDME takes V/kappa (1 + h*/h) = 1.5 V/kappa
        // at h = 2 h*, the corrected one V/kappa. The 11^3 zero-flux grid
        // differs from the periodic, N -> infinity derivation by a few per
        // cent (measured with 4000 trajectories: 0.966 +- 0.015 and
        // 1.458 +- 0.023), hence the 5 % allowance on top of the noise.
        let runs = 2000;
        let (fixed, sem) = association_time(true, runs);
        assert!(
            (fixed - 1.0).abs() < 0.05 + 3.0 * sem,
            "corrected: {fixed} +- {sem}"
        );
        let (plain, sem_plain) = association_time(false, runs);
        assert!(
            (plain - 1.5).abs() < 0.05 + 3.0 * sem_plain,
            "uncorrected: {plain} +- {sem_plain}"
        );
        assert!(
            plain - fixed > 0.05 + 3.0 * (sem + sem_plain),
            "{fixed} vs {plain}"
        );
    }
}
