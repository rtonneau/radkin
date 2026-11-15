//! Multi-grid RDME: every species on its own level of the voxel hierarchy
//! (ADR 0005; Cao and Erban 2014, Erban and Winkelmann 2024).
//!
//! A species of species level `L` holds its counts on the voxels of level
//! `L` (2^L x 2^L x 2^L blocks of fine voxels) and hops there, with
//! propensity `D / h_L^2` per molecule and face. A reaction belongs to the
//! level of its coarsest reactant and fires in a voxel of that level, with
//! the propensity of that level's model; the counts of a finer reactant are
//! summed over the fine voxels the coarse voxel contains (Hellander and
//! Hellander 2020, III.C.4: the pair reacts as if both occupied the coarser
//! voxel).
//!
//! Firing a reaction:
//!
//! - a consumed molecule of a finer reactant is taken from one of the
//!   contained voxels, in proportion to their counts;
//! - a product goes where the finest consumed molecule was: to the voxel of
//!   its own level that contains that spot, or, if the product lives on a
//!   finer level than that molecule, uniformly to one of the voxels of its
//!   level inside it (Erban and Winkelmann 2024, eq. 20). Without a consumed
//!   molecule the reaction voxel is the spot.
//!
//! The follower of a reversible pair (ADR 0006) carries the mesh factor of
//! the level where its bimolecular direction fires, so K_eq is kept when the
//! two directions belong to different levels.
//!
//! Each species also keeps its counts summed over the voxels of every
//! coarser level whose reactions read it, so a propensity reads its counts
//! directly. The event queue covers the voxels of every level
//! (next-subvolume method over all levels).

use crate::error::{RadkinError, Result};
use crate::ssa::model::{Law, Slot, SsaModel};
use crate::ssa::rdme::{EventQueue, RdmeSetup, ReversiblePair, face_weight, scale_down, strides};
use crate::ssa::rng::SsaRng;

/// The levels of species and reactions at one global coarsening level.
#[derive(Debug, Clone)]
pub(crate) struct Layout {
    /// Level every species lives on: `max(global level, species level)`.
    pub(crate) level: Vec<usize>,
    /// Highest level whose reactions read each species (at least its own).
    pub(crate) top: Vec<usize>,
    /// Level of every reaction: the coarsest level of its reactants.
    pub(crate) react_level: Vec<usize>,
    /// Reactions of every level.
    pub(crate) by_level: Vec<Vec<usize>>,
    /// `readers[l][s]`: reactions of level `l` whose propensity reads `s`.
    pub(crate) readers: Vec<Vec<Vec<usize>>>,
    /// Species read by the reactions of every level.
    pub(crate) read_at: Vec<Vec<usize>>,
    /// Diffusing species of every level.
    pub(crate) hoppers: Vec<Vec<usize>>,
    /// Propensity model of every level: the base model of that level, except
    /// that the follower of a shared reversible pair carries the mesh factor
    /// of the level where its bimolecular direction fires (ADR 0006).
    pub(crate) models: Vec<SsaModel>,
}

/// The levels of species and reactions, shared by the trajectories.
#[derive(Debug, Clone)]
pub struct MultigridSetup {
    pub(crate) base: RdmeSetup,
    /// Layout at every global coarsening level, `0..=max_global`; the
    /// layout of level 0 holds the species levels.
    pub(crate) layouts: Vec<Layout>,
    /// First site of every level, then the number of sites.
    pub(crate) offset: Vec<usize>,
    /// `parent[l][v]`: the voxel of level `l + 1` that contains voxel `v` of level `l`.
    parent: Vec<Vec<u32>>,
}

/// Dynamic species read by the propensity of reaction `r`.
fn reads(r: &crate::ssa::model::Reaction) -> Vec<usize> {
    let mut out: Vec<usize> = r.need.iter().map(|(s, _)| *s).collect();
    match &r.law {
        Law::MassAction { reactants, .. } => out.extend(reactants.iter().map(|(s, _)| *s)),
        Law::Michaelis {
            enzyme, substrate, ..
        } => {
            for s in [enzyme, substrate] {
                if let Slot::Dyn(i) = s {
                    out.push(*i);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

impl Layout {
    fn new(base: &RdmeSetup, level: Vec<usize>) -> Self {
        let n_sp = base.n_sp;
        let depth = base.max_levels() as usize;
        let model = base.model();
        let react_level: Vec<usize> = model
            .reactions
            .iter()
            .map(|r| {
                let read = reads(r);
                if read.is_empty() {
                    // No dynamic reactant: the finest level among the products.
                    r.delta.iter().map(|(s, _)| level[*s]).min().unwrap_or(0)
                } else {
                    read.iter().map(|s| level[*s]).max().unwrap_or(0)
                }
            })
            .collect();
        let mut by_level = vec![Vec::new(); depth + 1];
        let mut readers = vec![vec![Vec::new(); n_sp]; depth + 1];
        let mut top = level.clone();
        for (j, r) in model.reactions.iter().enumerate() {
            let l = react_level[j];
            by_level[l].push(j);
            for s in reads(r) {
                readers[l][s].push(j);
                top[s] = top[s].max(l);
            }
        }
        let read_at = (0..=depth)
            .map(|l| (0..n_sp).filter(|s| !readers[l][*s].is_empty()).collect())
            .collect();
        let hoppers = (0..=depth)
            .map(|l| {
                (0..n_sp)
                    .filter(|s| level[*s] == l && base.diffusion[*s] > 0.0)
                    .collect()
            })
            .collect();
        let mut models: Vec<SsaModel> = base.levels.iter().map(|lv| lv.model.clone()).collect();
        for p in base.reversible_pairs() {
            let ReversiblePair::Shared {
                bimolecular,
                follower,
                ..
            } = *p
            else {
                continue;
            };
            let (lb, lf) = (react_level[bimolecular], react_level[follower]);
            if lb != lf {
                // The follower fires on level `lf` with the factor of that
                // level: swap it for the factor of level `lb`.
                let by = base.pole_distance(lb, bimolecular) / base.pole_distance(lf, bimolecular);
                scale_down(&mut models[lf], follower, by);
            }
        }
        Self {
            level,
            top,
            react_level,
            by_level,
            readers,
            read_at,
            hoppers,
            models,
        }
    }
}

impl MultigridSetup {
    /// `base` with the species levels `levels` (state order), and global
    /// coarsening allowed up to level `max_global`; every level must exist
    /// in `base`.
    pub fn new(base: RdmeSetup, levels: &[u32], max_global: u32) -> Result<Self> {
        let n_sp = base.n_sp;
        let depth = base.max_levels() as usize;
        if levels.len() != n_sp {
            return Err(RadkinError::Value(format!(
                "rdme: {} species levels for {n_sp} dynamic species",
                levels.len()
            )));
        }
        if let Some(l) = levels
            .iter()
            .chain(std::iter::once(&max_global))
            .find(|l| **l as usize > depth)
        {
            return Err(RadkinError::Value(format!(
                "rdme: level {l} above the {depth} levels of the grid"
            )));
        }
        let species_level: Vec<usize> = levels.iter().map(|l| *l as usize).collect();
        let layouts = (0..=max_global as usize)
            .map(|g| Layout::new(&base, species_level.iter().map(|l| (*l).max(g)).collect()))
            .collect();
        let mut offset = vec![0];
        for lv in &base.levels {
            offset.push(offset.last().unwrap() + lv.grid.n_voxels());
        }
        let parent = (0..depth)
            .map(|l| {
                let (fine, coarse) = (&base.levels[l].grid, &base.levels[l + 1].grid);
                (0..fine.n_voxels())
                    .map(|v| {
                        let c = fine.coords(v);
                        coarse.index(c[0] / 2, c[1] / 2, c[2] / 2) as u32
                    })
                    .collect()
            })
            .collect();
        Ok(Self {
            base,
            layouts,
            offset,
            parent,
        })
    }

    /// Global coarsening levels allowed above the fine grid.
    pub fn max_global(&self) -> u32 {
        self.layouts.len() as u32 - 1
    }

    /// The single-level setup it extends.
    pub fn base(&self) -> &RdmeSetup {
        &self.base
    }

    /// The first reaction whose reactants and products are not all on the
    /// level of the reaction.
    pub fn cross_level_reaction(&self) -> Option<usize> {
        let model = self.base.model();
        let lay = &self.layouts[0];
        (0..model.reactions.len()).find(|j| {
            let l = lay.react_level[*j];
            let r = &model.reactions[*j];
            reads(r)
                .iter()
                .chain(r.delta.iter().map(|(s, _)| s))
                .any(|s| lay.level[*s] != l)
        })
    }

    /// Level and voxel of site `i`.
    pub(crate) fn site_of(&self, i: usize) -> (usize, usize) {
        let l = self.offset.partition_point(|o| *o <= i) - 1;
        (l, i - self.offset[l])
    }

    /// Voxel of level `to` that contains voxel `v` of level `from <= to`.
    pub(crate) fn ancestor(&self, mut v: usize, from: usize, to: usize) -> usize {
        for parent in &self.parent[from..to] {
            v = parent[v] as usize;
        }
        v
    }
}

/// One trajectory of the multi-grid RDME.
#[derive(Debug, Clone)]
pub struct Multigrid<'a> {
    pub(crate) setup: &'a MultigridSetup,
    pub t: f64,
    /// `agg[s][k][v]`: molecules of `s` in voxel `v` of level `level[s] + k`;
    /// `k = 0` holds its own counts, `k` goes up to `top[s] - level[s]`.
    pub(crate) agg: Vec<Vec<Vec<i64>>>,
    /// Total reaction propensity of every site.
    react: Vec<f64>,
    /// `sum_s D_s n_s` over the species living on every site.
    dsum: Vec<f64>,
    queue: EventQueue,
    pub(crate) rng: SsaRng,
    /// Scratch: counts seen by the reactions of one site, species-indexed.
    pub(crate) buf: Vec<i64>,
    /// Scratch: propensities of the reactions of one site.
    a: Vec<f64>,
    /// Scratch: sites whose rates changed.
    touched: Vec<usize>,
    pub n_reactions: u64,
    pub n_hops: u64,
    /// Events selected on a site whose rate had become zero (round-off of
    /// the incremental rate updates); they change nothing.
    pub n_null: u64,
    /// Current global coarsening level: every species lives on at least it.
    g: usize,
    /// Threshold of the dispersion test; `None`: global coarsening off.
    z_max: Option<f64>,
    /// `(t, global level)` after every change of global level.
    pub level_history: Vec<(f64, u32)>,
    /// Events fired since the last coarsening test.
    events_since_test: u64,
    /// Tau-leaping state; `None`: exact.
    pub(crate) tau: Option<Box<crate::ssa::rdme_multigrid_tau::MultigridTau>>,
    /// The event queue does not match the counts (set by a leap).
    pub(crate) queue_stale: bool,
}

/// Empty count arrays of every species for `lay`: its own level, then the
/// coarser levels up to `top`.
fn empty_agg(setup: &MultigridSetup, lay: &Layout) -> Vec<Vec<Vec<i64>>> {
    (0..setup.base.n_sp)
        .map(|s| {
            (lay.level[s]..=lay.top[s])
                .map(|l| vec![0; setup.base.levels[l].grid.n_voxels()])
                .collect()
        })
        .collect()
}

impl<'a> Multigrid<'a> {
    /// An empty grid at time `t`.
    pub fn new(setup: &'a MultigridSetup, t: f64, rng: SsaRng) -> Self {
        let base = &setup.base;
        let agg = empty_agg(setup, &setup.layouts[0]);
        let max_by_level = setup
            .layouts
            .iter()
            .flat_map(|lay| lay.by_level.iter().map(Vec::len))
            .max()
            .unwrap_or(0);
        let mut r = Self {
            setup,
            g: 0,
            z_max: None,
            level_history: Vec::new(),
            events_since_test: 0,
            tau: None,
            queue_stale: false,
            t,
            agg,
            react: Vec::new(),
            dsum: Vec::new(),
            queue: EventQueue::default(),
            rng,
            buf: vec![0; base.n_sp],
            a: vec![0.0; max_by_level],
            touched: Vec::new(),
            n_reactions: 0,
            n_hops: 0,
            n_null: 0,
        };
        r.rebuild();
        r
    }

    /// Reaction and hop events fired.
    pub fn n_events(&self) -> u64 {
        self.n_reactions + self.n_hops
    }

    /// Layout of the current global level.
    pub(crate) fn lay(&self) -> &'a Layout {
        &self.setup.layouts[self.g]
    }

    /// Current global coarsening level.
    pub fn global_level(&self) -> u32 {
        self.g as u32
    }

    /// Level species `s` currently lives on.
    pub fn level_of(&self, s: usize) -> u32 {
        self.lay().level[s] as u32
    }

    /// Turn global coarsening on with the dispersion threshold `z` (`None`: off).
    pub fn set_coarsening(&mut self, z: Option<f64>) {
        self.z_max = z;
    }

    /// Move to global level `g`: the species whose level changes are merged
    /// (counts summed) or split (each molecule to one of the 8 children at
    /// random, the uniform multinomial split), one level at a time.
    fn set_global(&mut self, g: usize) {
        let setup = self.setup;
        let (old, new) = (self.lay(), &setup.layouts[g]);
        let mut own: Vec<Vec<i64>> = self
            .agg
            .iter_mut()
            .map(|a| std::mem::take(&mut a[0]))
            .collect();
        for (s, counts) in own.iter_mut().enumerate() {
            let (mut l, to) = (old.level[s], new.level[s]);
            while l < to {
                let mut up = vec![0; setup.base.levels[l + 1].grid.n_voxels()];
                for (v, n) in counts.iter().enumerate() {
                    up[setup.parent[l][v] as usize] += n;
                }
                *counts = up;
                l += 1;
            }
            while l > to {
                let fine = &setup.base.levels[l - 1].grid;
                let coarse = &setup.base.levels[l].grid;
                let mut down = vec![0; fine.n_voxels()];
                for (v, n) in counts.iter().enumerate() {
                    let c = coarse.coords(v);
                    for _ in 0..*n {
                        let k = self.rng.below(8) as usize;
                        down[fine.index(
                            2 * c[0] + (k & 1),
                            2 * c[1] + ((k >> 1) & 1),
                            2 * c[2] + ((k >> 2) & 1),
                        )] += 1;
                    }
                }
                *counts = down;
                l -= 1;
            }
        }
        self.agg = empty_agg(setup, new);
        for (a, counts) in self.agg.iter_mut().zip(own) {
            a[0] = counts;
        }
        self.g = g;
        self.level_history.push((self.t, g as u32));
        self.rebuild();
    }

    /// Whether the species that would merge at the next global level (those
    /// living on the current one) are consistent with a uniform split over
    /// every 2x2x2 block (the test of [`crate::ssa::rdme_coarsen`]). False
    /// at the last global level.
    pub fn is_uniform_in_blocks(&self, z: f64) -> bool {
        let setup = self.setup;
        let g = self.g;
        if g + 1 >= setup.layouts.len() {
            return false;
        }
        let lay = self.lay();
        let species: Vec<usize> = (0..setup.base.n_sp)
            .filter(|s| lay.level[*s] == g)
            .collect();
        crate::ssa::rdme_coarsen::blocks_are_uniform(
            &setup.base.levels[g].grid,
            &setup.base.levels[g + 1].grid,
            &species,
            |v, s| self.agg[s][0][v],
            z,
        )
    }

    /// Drop as many global levels as the dispersion test allows. Returns
    /// whether the level changed.
    pub fn coarsen_while_uniform(&mut self) -> bool {
        let Some(z) = self.z_max else {
            return false;
        };
        let before = self.g;
        while self.is_uniform_in_blocks(z) {
            self.set_global(self.g + 1);
        }
        self.g != before
    }

    /// Back to global level 0: every species on its own species level.
    pub fn refine_to_fine(&mut self) {
        while self.g > 0 {
            self.set_global(self.g - 1);
        }
    }

    /// Run the coarsening test once `8 x voxels of the current global level`
    /// events fired since the last one.
    pub(crate) fn count_events(&mut self, fired: u64) {
        if self.z_max.is_none() {
            return;
        }
        self.events_since_test += fired;
        let interval = 8 * self.setup.base.levels[self.g].grid.n_voxels() as u64;
        if self.events_since_test >= interval {
            self.events_since_test = 0;
            self.coarsen_while_uniform();
        }
    }

    /// Sum the own counts of every species into its coarser levels.
    pub(crate) fn rebuild_aggregates(&mut self) {
        let setup = self.setup;
        let lay = self.lay();
        for (s, levels) in self.agg.iter_mut().enumerate() {
            let l0 = lay.level[s];
            for k in 1..levels.len() {
                let (lower, upper) = levels.split_at_mut(k);
                let (fine, coarse) = (&lower[k - 1], &mut upper[0]);
                coarse.iter_mut().for_each(|x| *x = 0);
                for (v, n) in fine.iter().enumerate() {
                    coarse[setup.parent[l0 + k - 1][v] as usize] += n;
                }
            }
        }
    }

    /// Fill `buf` with the counts the reactions of voxel `v` of level `l` read.
    pub(crate) fn load(&mut self, l: usize, v: usize) {
        let lay = self.lay();
        for &s in &lay.read_at[l] {
            self.buf[s] = self.agg[s][l - lay.level[s]][v];
        }
    }

    /// Total reaction propensity of voxel `v` of level `l`.
    fn react_of(&mut self, l: usize, v: usize) -> f64 {
        self.load(l, v);
        let model = &self.lay().models[l];
        self.lay().by_level[l]
            .iter()
            .map(|j| model.propensity(*j, &self.buf))
            .sum()
    }

    /// Propensity of the reactions of level `l` that read `s`, in voxel `v`.
    fn readers_sum(&mut self, l: usize, v: usize, s: usize) -> f64 {
        self.load(l, v);
        let model = &self.lay().models[l];
        self.lay().readers[l][s]
            .iter()
            .map(|j| model.propensity(*j, &self.buf))
            .sum()
    }

    /// `sum_s D_s n_s` of the species living on voxel `v` of level `l`.
    fn dsum_of(&self, l: usize, v: usize) -> f64 {
        let setup = self.setup;
        let lay = self.lay();
        lay.hoppers[l]
            .iter()
            .map(|s| self.agg[*s][0][v] as f64 * setup.base.diffusion[*s])
            .sum()
    }

    /// Total event rate of site `i` (level `l`, voxel `v`).
    fn rate(&self, i: usize, l: usize, v: usize) -> f64 {
        let lv = &self.setup.base.levels[l];
        self.react[i] + self.dsum[i] * face_weight(&lv.grid, &lv.inv_h2, v)
    }

    /// Recompute every rate and draw every next event time.
    pub(crate) fn rebuild(&mut self) {
        self.rebuild_aggregates();
        let setup = self.setup;
        let n_sites = *setup.offset.last().unwrap();
        self.react = vec![0.0; n_sites];
        self.dsum = vec![0.0; n_sites];
        let mut keys = Vec::with_capacity(n_sites);
        for i in 0..n_sites {
            let (l, v) = setup.site_of(i);
            self.react[i] = self.react_of(l, v);
            self.dsum[i] = self.dsum_of(l, v);
            let rate = self.rate(i, l, v);
            keys.push(if rate > 0.0 {
                self.t + self.rng.exp1() / rate
            } else {
                f64::INFINITY
            });
        }
        self.queue = EventQueue::build(keys);
        self.queue_stale = false;
    }

    /// Set every count to a Poisson draw of mean `c N_A V` for the uniform
    /// concentrations `y` [µmol/L], `V` the voxel volume of the species' level.
    pub fn fill_poisson(&mut self, y: &[f64]) {
        let setup = self.setup;
        let lay = self.lay();
        assert_eq!(y.len(), setup.base.n_sp, "rdme: state length");
        for (s, c) in y.iter().enumerate() {
            let k = setup.base.levels[lay.level[s]].model.molecules_per_umol();
            for n in self.agg[s][0].iter_mut() {
                *n = self.rng.poisson(c * k);
            }
        }
        self.rebuild();
    }

    /// Add `(fine voxel, species, molecules)` entries, each to the voxel of
    /// its species' level that contains the fine voxel. A globally coarsened
    /// grid is refined to global level 0 first.
    pub fn inject(&mut self, molecules: &[(usize, usize, i64)]) {
        self.refine_to_fine();
        self.events_since_test = 0;
        let setup = self.setup;
        let lay = self.lay();
        for &(v, s, n) in molecules {
            let w = setup.ancestor(v, 0, lay.level[s]);
            self.agg[s][0][w] += n;
        }
        self.rebuild();
    }

    /// Molecules of every species in the domain.
    pub fn totals(&self) -> Vec<i64> {
        self.agg
            .iter()
            .map(|levels| levels[0].iter().sum())
            .collect()
    }

    /// Counts of species `s` on the voxels of its own level.
    pub fn counts(&self, s: usize) -> &[i64] {
        &self.agg[s][0]
    }

    /// The field of every species on the fine grid, species-major
    /// (`out[s * n_fine + v]`), in units of `1 / 8^max_levels` molecule: a
    /// voxel of level `L` spreads its count evenly over its `8^L` fine voxels.
    pub fn fine_field(&self) -> Vec<i64> {
        let setup = self.setup;
        let lay = self.lay();
        let n_fine = setup.base.grid().n_voxels();
        let max = setup.base.max_levels() as usize;
        let mut out = vec![0; setup.base.n_sp * n_fine];
        for (s, levels) in self.agg.iter().enumerate() {
            let l = lay.level[s];
            let scale = 8i64.pow((max - l) as u32);
            for v in 0..n_fine {
                out[s * n_fine + v] = levels[0][setup.ancestor(v, 0, l)] * scale;
            }
        }
        out
    }

    /// Add `d` molecules of `s` to voxel `w` of its level, keep the sums of
    /// the coarser levels, and mark the sites whose rates it changes.
    fn change(&mut self, s: usize, w: usize, d: i64) {
        let setup = self.setup;
        let lay = self.lay();
        let l0 = lay.level[s];
        let mut v = w;
        for k in 0..self.agg[s].len() {
            let l = l0 + k;
            if k > 0 {
                v = setup.parent[l - 1][v] as usize;
            }
            self.agg[s][k][v] += d;
            debug_assert!(self.agg[s][k][v] >= 0);
            if k == 0 || !lay.readers[l][s].is_empty() {
                self.touched.push(setup.offset[l] + v);
            }
        }
    }

    /// A voxel of level `level[s]` inside voxel `v` of level `l`, picked in
    /// proportion to the molecules of `s` it holds.
    pub(crate) fn pick_inside(&mut self, s: usize, l: usize, mut v: usize) -> usize {
        let setup = self.setup;
        let lay = self.lay();
        let l0 = lay.level[s];
        for cur in (l0 + 1..=l).rev() {
            let fine = &setup.base.levels[cur - 1].grid;
            let c = setup.base.levels[cur].grid.coords(v);
            let mut x = self.rng.below(self.agg[s][cur - l0][v] as u64) as i64;
            let mut chosen = None;
            'kids: for n in 0..8 {
                let child = fine.index(
                    2 * c[0] + (n & 1),
                    2 * c[1] + ((n >> 1) & 1),
                    2 * c[2] + ((n >> 2) & 1),
                );
                let m = self.agg[s][cur - 1 - l0][child];
                if x < m {
                    chosen = Some(child);
                    break 'kids;
                }
                x -= m;
            }
            v = chosen.expect("rdme: the sums of the levels disagree");
        }
        v
    }

    /// A voxel of level `to` inside voxel `v` of level `from > to`, uniformly.
    pub(crate) fn uniform_inside(&mut self, v: usize, from: usize, to: usize) -> usize {
        let levels = &self.setup.base.levels;
        let c = levels[from].grid.coords(v);
        let scale = 1usize << (from - to);
        let mut p = [0; 3];
        for a in 0..3 {
            p[a] = c[a] * scale + self.rng.below(scale as u64) as usize;
        }
        levels[to].grid.index(p[0], p[1], p[2])
    }

    /// Fire a reaction of voxel `v` of level `l`; `frac` in `[0, 1)` selects it.
    fn fire_reaction(&mut self, l: usize, v: usize, frac: f64) {
        let setup = self.setup;
        let lay = self.lay();
        let model = &lay.models[l];
        let site = setup.offset[l] + v;
        self.load(l, v);
        let reactions = &lay.by_level[l];
        let mut total = 0.0;
        for (a, j) in self.a.iter_mut().zip(reactions) {
            *a = model.propensity(*j, &self.buf);
            total += *a;
        }
        if total <= 0.0 {
            self.react[site] = 0.0;
            self.n_null += 1;
            return;
        }
        let target = frac * total;
        let (mut acc, mut chosen) = (0.0, reactions[0]);
        for (a, j) in self.a.iter().zip(reactions) {
            if *a > 0.0 {
                acc += a;
                chosen = *j;
                if target < acc {
                    break;
                }
            }
        }
        self.touched.clear();
        // The spot of the finest consumed molecule: (level, voxel).
        let mut spot = (l, v);
        let delta = &model.reactions[chosen].delta;
        for &(s, d) in delta.iter().filter(|(_, d)| *d < 0) {
            for _ in 0..-d {
                let w = self.pick_inside(s, l, v);
                self.change(s, w, -1);
                if lay.level[s] < spot.0 {
                    spot = (lay.level[s], w);
                }
            }
        }
        for &(s, d) in delta.iter().filter(|(_, d)| *d > 0) {
            let ls = lay.level[s];
            for _ in 0..d {
                let w = if ls >= spot.0 {
                    setup.ancestor(spot.1, spot.0, ls)
                } else {
                    self.uniform_inside(spot.1, spot.0, ls)
                };
                self.change(s, w, 1);
            }
        }
        self.touched.sort_unstable();
        self.touched.dedup();
        let touched = std::mem::take(&mut self.touched);
        for &i in &touched {
            let (li, vi) = setup.site_of(i);
            self.react[i] = self.react_of(li, vi);
            self.dsum[i] = self.dsum_of(li, vi);
            if i != site {
                self.reschedule(i, li, vi);
            }
        }
        self.touched = touched;
        if !self.touched.contains(&site) {
            self.react[site] = self.react_of(l, v);
        }
        self.n_reactions += 1;
    }

    /// Fire a hop out of voxel `v` of level `l`; `x` in `[0, dsum)` selects
    /// the species.
    fn fire_hop(&mut self, l: usize, v: usize, x: f64) {
        let setup = self.setup;
        let lay = self.lay();
        let site = setup.offset[l] + v;
        let (mut acc, mut species) = (0.0, None);
        for &s in &lay.hoppers[l] {
            let n = self.agg[s][0][v];
            if n > 0 {
                acc += n as f64 * setup.base.diffusion[s];
                species = Some(s);
                if x < acc {
                    break;
                }
            }
        }
        let Some(s) = species else {
            self.dsum[site] = 0.0;
            self.n_null += 1;
            return;
        };
        let lv = &setup.base.levels[l];
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
            self.n_null += 1;
            return;
        }
        // Move the molecule on its level and on every coarser level it
        // leaves a voxel of; update the reactions that read it there.
        let (mut a, mut b) = (v, dest);
        for k in 0..self.agg[s].len() {
            let lk = l + k;
            if k > 0 {
                a = setup.parent[lk - 1][a] as usize;
                b = setup.parent[lk - 1][b] as usize;
                if a == b {
                    break;
                }
            }
            let has_readers = !lay.readers[lk][s].is_empty();
            let before = if has_readers {
                (self.readers_sum(lk, a, s), self.readers_sum(lk, b, s))
            } else {
                (0.0, 0.0)
            };
            self.agg[s][k][a] -= 1;
            self.agg[s][k][b] += 1;
            let (ia, ib) = (setup.offset[lk] + a, setup.offset[lk] + b);
            if has_readers {
                let after = (self.readers_sum(lk, a, s), self.readers_sum(lk, b, s));
                self.react[ia] = (self.react[ia] + (after.0 - before.0)).max(0.0);
                self.react[ib] = (self.react[ib] + (after.1 - before.1)).max(0.0);
            }
            if k == 0 {
                self.dsum[ia] = self.dsum_of(lk, a);
                self.dsum[ib] = self.dsum_of(lk, b);
                self.reschedule(ib, lk, b);
            } else if has_readers {
                self.reschedule(ia, lk, a);
                self.reschedule(ib, lk, b);
            }
        }
        self.n_hops += 1;
    }

    /// New next-event time of site `i` from its current rates.
    fn reschedule(&mut self, i: usize, l: usize, v: usize) {
        let rate = self.rate(i, l, v);
        let key = if rate > 0.0 {
            self.t + self.rng.exp1() / rate
        } else {
            f64::INFINITY
        };
        self.queue.set(i, key);
    }

    /// Fire the next event if it happens before `t_stop`. Returns whether
    /// one did.
    pub fn step(&mut self, t_stop: f64) -> bool {
        let (i, t_next) = self.queue.top();
        if t_next >= t_stop {
            return false;
        }
        self.t = t_next;
        let (l, v) = self.setup.site_of(i);
        let lv = &self.setup.base.levels[l];
        let hop = self.dsum[i] * face_weight(&lv.grid, &lv.inv_h2, v);
        let u = self.rng.uniform() * (self.react[i] + hop);
        if u < hop {
            self.fire_hop(l, v, u / hop * self.dsum[i]);
        } else {
            let frac = ((u - hop) / self.react[i]).clamp(0.0, 1.0 - f64::EPSILON);
            self.fire_reaction(l, v, frac);
        }
        self.reschedule(i, l, v);
        true
    }

    /// Panics unless the sums of the coarser levels and the rates of every
    /// site match a recomputation from the own counts.
    #[cfg(test)]
    fn check_consistency(&mut self) {
        let before = self.agg.clone();
        self.rebuild_aggregates();
        assert_eq!(before, self.agg, "level sums");
        let setup = self.setup;
        for i in 0..self.react.len() {
            let (l, v) = setup.site_of(i);
            let exact = self.react_of(l, v);
            assert!(
                (self.react[i] - exact).abs() <= 1e-9 * exact.max(1.0),
                "site {i}: {} vs {exact}",
                self.react[i]
            );
            assert_eq!(self.dsum[i], self.dsum_of(l, v), "site {i}");
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
    use crate::matrices::SimMatrices;
    use crate::ssa::rdme::Rdme;
    use crate::voxel::grid::VoxelGrid;

    fn setup(
        toml: &str,
        dims: [usize; 3],
        h: f64,
        d: f64,
        depth: u32,
        levels: &[(&str, u32)],
    ) -> (SimEnv, MultigridSetup) {
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new(dims, [h; 3], [0.0; 3]).unwrap();
        let base = RdmeSetup::new(&mat, &grid, vec![d; mat.n_dyn], depth, false).unwrap();
        let mut l = vec![0; mat.n_dyn];
        for (name, level) in levels {
            l[env.species.index_of(name).unwrap()] = *level;
        }
        (env, MultigridSetup::new(base, &l, 0).unwrap())
    }

    const DECAY: &str = r#"
        [[reactions.k_reaction]]
        reaction = "H_r -> H2"
        k_value = 2.0e6
    "#;

    #[test]
    fn levels_of_reactions_follow_the_coarsest_reactant() {
        let toml = r#"
            [[reactions.k_reaction]]
            reaction = "e_aq + OH_r -> OH_minus"
            k_value = 3.0e10
            [[reactions.k_reaction]]
            reaction = "H_r -> H2"
            k_value = 2.0e6
        "#;
        let (_, s) = setup(toml, [4, 4, 4], 10e-9, 5e-9, 2, &[("OH_r", 2)]);
        assert_eq!(s.layouts[0].react_level, [2, 0]);
        assert_eq!(s.cross_level_reaction(), Some(0));
        let (_, same) = setup(toml, [4, 4, 4], 10e-9, 5e-9, 2, &[]);
        assert_eq!(same.cross_level_reaction(), None);
        let (_, decay) = setup(DECAY, [4, 4, 4], 10e-9, 5e-9, 2, &[("H_r", 1), ("H2", 1)]);
        assert_eq!(decay.cross_level_reaction(), None);
        assert_eq!(decay.layouts[0].react_level, [1]);
    }

    #[test]
    fn a_reversible_pair_takes_the_factor_of_its_bimolecular_level() {
        // H2O <-> H3O+ + OH-: the forward direction fires on the finest
        // product level, the backward one on the coarsest reactant level.
        // With OH- on level 1, the forward direction on level 0 must carry
        // the factor of level 1, where the backward one fires (ADR 0006).
        let toml = r#"
            [[reactions.equilibrium]]
            reaction = "H2O <-> H3O_plus + OH_minus"
            k_forward = 2.12e-5
            k_backward = 1.18e11
        "#;
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        // h* = 4.95 nm (D = 5e-9 each): factors 1.98 (10 nm) and 1.33 (20 nm).
        let h = 10e-9;
        let grid = VoxelGrid::new([4, 4, 4], [h; 3], [0.0; 3]).unwrap();
        let oh = env.species.index_of("OH_minus").unwrap();
        let mut levels = vec![0; mat.n_dyn];
        levels[oh] = 1;
        let make = |correct: bool, levels: &[u32], global: u32| {
            let base = RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], 1, correct).unwrap();
            MultigridSetup::new(base, levels, global).unwrap()
        };
        let c_of = |m: &SsaModel, j: usize| match &m.reactions[j].law {
            Law::MassAction { c, .. } => *c,
            other => panic!("{other:?}"),
        };
        // `c` of reaction `j` on level `l` of the layout of global level `g`.
        let c = |s: &MultigridSetup, g: usize, l: usize, j: usize| c_of(&s.layouts[g].models[l], j);
        let (plain, fixed) = (make(false, &levels, 1), make(true, &levels, 1));
        assert_eq!(fixed.layouts[0].react_level, [0, 1]);
        let hs = fixed.base.critical_spacing(1).unwrap();
        let f1 = 1.0 / (1.0 - hs / (2.0 * h));
        for (j, l) in [(0, 0), (1, 1)] {
            let got = c(&fixed, 0, l, j) / c(&plain, 0, l, j);
            assert!((got - f1).abs() < 1e-12 * f1, "reaction {j}: {got} vs {f1}");
        }
        // The base model of level 0 keeps the factor of level 0.
        let f0 = 1.0 / (1.0 - hs / h);
        let base0 = c_of(&fixed.base.levels[0].model, 0);
        assert!((base0 / c(&plain, 0, 0, 0) - f0).abs() < 1e-12 * f0);
        assert!(f0 > 1.4 * f1);
        // At global level 1 both directions fire on level 1: nothing to move.
        assert_eq!(fixed.layouts[1].react_level, [1, 1]);
        assert_eq!(fixed.layouts[1].models[1], fixed.base.levels[1].model);
        // Same levels: the layout models are the base ones.
        let same = make(true, &vec![0; mat.n_dyn], 0);
        assert_eq!(same.layouts[0].models[0], same.base.levels[0].model);
        assert_eq!(same.layouts[0].models[1], same.base.levels[1].model);
    }

    #[test]
    fn a_coarse_species_diffuses_on_its_level() {
        // 8x2x2 fine voxels of 50 nm; level 1 is a row of 4 voxels of
        // 100 nm with faces along x only. The mean count of a hop process
        // started in voxel 0 solves dm/dt = (D/h^2) L m exactly:
        // m_j = N/n [1 + 2 sum_k cos(pi k/2n) cos(pi k (j+1/2)/n) e^(-lambda_k t)],
        // lambda_k = 2 (1 - cos(pi k/n)) D/h^2.
        let (env, s) = setup(
            "[reactions.radiolytic]\ne_aq = 2.6",
            [8, 2, 2],
            50e-9,
            5e-9,
            1,
            &[("e_aq", 1)],
        );
        let e = env.species.index_of("e_aq").unwrap();
        let (n0, runs, t) = (400i64, 200u64, 1.0e-6);
        let rate = 5e-9 / (100e-9f64).powi(2);
        let n = 4.0;
        let expected: Vec<f64> = (0..4)
            .map(|j| {
                let mut m = 1.0;
                for k in 1..4 {
                    let k = k as f64;
                    let lam = 2.0 * (1.0 - (std::f64::consts::PI * k / n).cos()) * rate;
                    m += 2.0
                        * (std::f64::consts::PI * k * 0.5 / n).cos()
                        * (std::f64::consts::PI * k * (j as f64 + 0.5) / n).cos()
                        * (-lam * t).exp();
                }
                n0 as f64 / n * m
            })
            .collect();
        let mut sum = [0.0; 4];
        for stream in 0..runs {
            let mut r = Multigrid::new(&s, 0.0, SsaRng::new(9, stream));
            // Fine voxel 0 lies in coarse voxel 0.
            r.inject(&[(0, e, n0)]);
            r.advance_to(t);
            assert_eq!(r.totals()[e], n0);
            assert_eq!(r.n_reactions, 0);
            for (j, x) in r.counts(e).iter().enumerate() {
                sum[j] += *x as f64;
            }
        }
        for j in 0..4 {
            let p = expected[j] / n0 as f64;
            let sem = (n0 as f64 * p * (1.0 - p) / runs as f64).sqrt();
            let mean = sum[j] / runs as f64;
            assert!(
                (mean - expected[j]).abs() <= 4.0 * sem,
                "voxel {j}: {mean} vs {}",
                expected[j]
            );
        }
    }

    #[test]
    fn a_coarse_reaction_has_the_macroscopic_rate() {
        // First-order decay on level 1: survival exp(-k t), as on the fine grid.
        let (env, s) = setup(DECAY, [4, 4, 4], 10e-9, 5e-9, 1, &[("H_r", 1), ("H2", 1)]);
        let (h, h2) = (
            env.species.index_of("H_r").unwrap(),
            env.species.index_of("H2").unwrap(),
        );
        let (n0, runs) = (500, 200);
        let mut sum = 0.0;
        for stream in 0..runs {
            let mut r = Multigrid::new(&s, 0.0, SsaRng::new(3, stream));
            r.inject(&[(21, h, n0)]);
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
    fn all_fine_levels_match_the_single_level_rdme() {
        // Same model, two algorithms: the decay survival agrees in law.
        let (env, s) = setup(DECAY, [3, 3, 3], 10e-9, 5e-9, 0, &[]);
        let h = env.species.index_of("H_r").unwrap();
        let (n0, runs) = (300, 200);
        let (mut a, mut b) = (0.0, 0.0);
        for stream in 0..runs {
            let mut m = Multigrid::new(&s, 0.0, SsaRng::new(5, stream));
            m.inject(&[(13, h, n0)]);
            m.advance_to(4e-7);
            a += m.totals()[h] as f64;
            let mut r = Rdme::new(s.base(), 0.0, SsaRng::new(6, stream));
            r.inject(&[(13, h, n0)]);
            r.advance_to(4e-7);
            b += r.totals()[h] as f64;
        }
        let p = (-2.0e6f64 * 4e-7).exp();
        let sem = (2.0 * n0 as f64 * p * (1.0 - p) / runs as f64).sqrt();
        assert!(((a - b) / runs as f64).abs() <= 4.0 * sem, "{a} vs {b}");
    }

    #[test]
    fn poisson_fill_and_fine_field_use_the_species_level() {
        let (env, s) = setup(DECAY, [4, 4, 4], 100e-9, 0.0, 2, &[("H_r", 2)]);
        let h = env.species.index_of("H_r").unwrap();
        let mut y = vec![0.0; s.base.n_sp];
        // 5 molecules per fine voxel: 320 per level-2 voxel.
        y[h] = 5.0 / s.base.model().molecules_per_umol();
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(2, 0));
        r.fill_poisson(&y);
        assert_eq!(r.counts(h).len(), 1);
        let total = r.totals()[h];
        assert!((total - 320).abs() < 4 * 18, "{total}");
        // Every fine voxel gets 1/64 of the level-2 count, in 1/8^2 units.
        let f = r.fine_field();
        let n_fine = 64;
        assert!(f[h * n_fine..(h + 1) * n_fine].iter().all(|x| *x == total));
        assert_eq!(
            f[h * n_fine..(h + 1) * n_fine].iter().sum::<i64>(),
            total * 64
        );
    }

    #[test]
    fn rates_and_level_sums_stay_consistent_across_levels() {
        let toml = r#"
            [[reactions.k_reaction]]
            reaction = "e_aq + OH_r -> OH_minus"
            k_value = 3.0e10
            [[reactions.k_reaction]]
            reaction = "OH_minus -> e_aq + OH_r"
            k_value = 1.0e6
            [[reactions.k_reaction]]
            reaction = "OH_r + OH_r -> H2O2"
            k_value = 5.5e9
        "#;
        let (env, s) = setup(toml, [4, 4, 4], 20e-9, 5e-9, 2, &[("OH_r", 1), ("H2O2", 2)]);
        let ix = |n: &str| env.species.index_of(n).unwrap();
        let (e, oh, ohm, h2o2) = (ix("e_aq"), ix("OH_r"), ix("OH_minus"), ix("H2O2"));
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(4, 0));
        r.inject(&[(0, e, 300), (63, oh, 300), (21, ohm, 50)]);
        for k in 1..=20 {
            r.advance_to(k as f64 * 2e-8);
            r.check_consistency();
            let t = r.totals();
            // e_aq + OH_minus and OH_r + OH_minus + 2 H2O2 are conserved.
            assert_eq!(t[e] + t[ohm], 350);
            assert_eq!(t[oh] + t[ohm] + 2 * t[h2o2], 350);
        }
        assert!(
            r.n_reactions > 100 && r.n_hops > 1000,
            "{} {}",
            r.n_reactions,
            r.n_hops
        );
        assert!(r.totals()[h2o2] > 0);
    }

    fn floor_setup(max_global: u32) -> (SimEnv, MultigridSetup) {
        let toml = r#"
            [[reactions.k_reaction]]
            reaction = "e_aq + OH_r -> OH_minus"
            k_value = 3.0e10
        "#;
        let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new([8, 4, 4], [20e-9; 3], [0.0; 3]).unwrap();
        let base = RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], 2, false).unwrap();
        let mut l = vec![0; mat.n_dyn];
        l[env.species.index_of("OH_r").unwrap()] = 1;
        (env, MultigridSetup::new(base, &l, max_global).unwrap())
    }

    #[test]
    fn global_levels_are_a_floor_and_conserve_every_species() {
        let (env, s) = floor_setup(2);
        let ix = |n: &str| env.species.index_of(n).unwrap();
        let (e, oh) = (ix("e_aq"), ix("OH_r"));
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(7, 0));
        let mut y = vec![0.0; s.base.n_sp];
        y[e] = 3.0 / s.base.model().molecules_per_umol();
        y[oh] = 2.0 / s.base.model().molecules_per_umol();
        r.fill_poisson(&y);
        let totals = r.totals();
        assert_eq!((r.level_of(e), r.level_of(oh)), (0, 1));
        r.set_global(1);
        assert_eq!((r.level_of(e), r.level_of(oh)), (1, 1));
        assert_eq!(r.counts(e).len(), 16);
        assert_eq!(r.totals(), totals);
        r.check_consistency();
        r.set_global(2);
        assert_eq!((r.level_of(e), r.level_of(oh)), (2, 2));
        assert_eq!(r.totals(), totals);
        r.check_consistency();
        // Refining stops at each species' own level.
        r.refine_to_fine();
        assert_eq!((r.level_of(e), r.level_of(oh)), (0, 1));
        assert_eq!(r.counts(e).len(), 128);
        assert_eq!(r.counts(oh).len(), 16);
        assert_eq!(r.totals(), totals);
        r.check_consistency();
        let history: Vec<u32> = r.level_history.iter().map(|(_, l)| *l).collect();
        assert_eq!(history, [1, 2, 1, 0]);
        // The fine field keeps the total in units of 1 / 8^2 molecule.
        let n_fine = 128;
        let f = r.fine_field();
        assert_eq!(
            f[e * n_fine..(e + 1) * n_fine].iter().sum::<i64>(),
            totals[e] * 64
        );
    }

    #[test]
    fn the_dispersion_test_reads_only_the_species_that_merge() {
        let (env, s) = floor_setup(2);
        let ix = |n: &str| env.species.index_of(n).unwrap();
        let (e, oh) = (ix("e_aq"), ix("OH_r"));
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(8, 0));
        // OH_r on level 1 is clustered in one voxel; at global level 0 only
        // e_aq (level 0) would merge, and it is absent.
        r.inject(&[(0, oh, 500)]);
        assert!(r.is_uniform_in_blocks(3.0));
        r.set_global(1);
        // Now OH_r merges into level 2: its cluster blocks the test.
        assert!(!r.is_uniform_in_blocks(3.0));
        r.refine_to_fine();
        r.inject(&[(5, e, 500)]);
        assert!(!r.is_uniform_in_blocks(3.0));
    }

    #[test]
    fn a_mixing_grid_coarsens_and_a_seed_refines_it() {
        let (env, s) = floor_setup(2);
        let ix = |n: &str| env.species.index_of(n).unwrap();
        let (e, oh, ohm) = (ix("e_aq"), ix("OH_r"), ix("OH_minus"));
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(9, 0));
        r.set_coarsening(Some(3.0));
        r.inject(&[(0, e, 2000), (127, oh, 2000)]);
        // h^2 / D = 80 ns: far from mixed after 2 ns.
        r.advance_to(2e-9);
        assert_eq!(r.global_level(), 0);
        r.advance_to(2e-5);
        assert_eq!(r.global_level(), 2, "history {:?}", r.level_history);
        r.check_consistency();
        let t = r.totals();
        assert_eq!(t[e] + t[ohm], 2000);
        assert_eq!(t[oh] + t[ohm], 2000);
        r.inject(&[(3, e, 10)]);
        assert_eq!(r.global_level(), 0);
        assert_eq!((r.level_of(e), r.level_of(oh)), (0, 1));
        assert_eq!(r.totals()[e], t[e] + 10);
        r.check_consistency();
    }

    #[test]
    fn injection_goes_to_the_containing_voxel() {
        let (env, s) = setup(DECAY, [4, 4, 4], 10e-9, 0.0, 1, &[("H_r", 1)]);
        let h = env.species.index_of("H_r").unwrap();
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(1, 0));
        // Fine (3, 3, 3) is in coarse (1, 1, 1), index 7.
        let fine = s.base.grid().index(3, 3, 3);
        r.inject(&[(fine, h, 5)]);
        assert_eq!(r.counts(h)[7], 5);
        assert_eq!(r.totals()[h], 5);
    }
}
