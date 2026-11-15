//! Coarsening of the RDME grid (global, one level at a time, as the
//! deterministic path: every 2x2x2 block of voxels becomes one voxel).
//!
//! The deterministic test (spread of the 8 children below a tolerance) can
//! never pass on molecule counts: Poisson noise alone exceeds it. The grid
//! drops a level when the counts of every block are statistically
//! consistent with a uniform distribution over its 8 voxels instead.
//!
//! For a block holding `S` molecules of a species, Pearson's statistic
//! `X = sum (c_i - S/8)^2 / (S/8)` has, for a uniform multinomial split,
//! the exact mean 7 and variance `14 (S - 1) / S`. Two conditions must hold
//! for every species:
//!
//! - aggregated over the blocks with `S >= 1`:
//!   `(sum X - 7 B) / sqrt(sum variance) <= z`;
//! - for each block with `S >= 16`:
//!   `(X - 7) / sqrt(variance) <= z + sqrt(2 ln B16)`, `B16` being the
//!   number of such blocks. This catches a localized cluster that the
//!   aggregate would dilute in a large uniform background.
//!
//! `z = 30 * coarsen_tol_factor` (3 at the default factor 0.1).
//!
//! Merging sums the 8 counts. Refining (before a seed injection) sends each
//! molecule of a voxel to one of its 8 children at random, which is the
//! uniform multinomial split.

use crate::ssa::rdme::Rdme;
use crate::voxel::grid::VoxelGrid;

/// `z` of the dispersion test for a `coarsen_tol_factor`.
pub fn dispersion_threshold(coarsen_tol_factor: f64) -> f64 {
    30.0 * coarsen_tol_factor
}

/// Blocks with at least this many molecules are also tested one by one.
const BLOCK_TEST_MIN: i64 = 16;

/// Flat indices of the 8 children of the coarse voxel `(i, j, k)` in `fine`.
fn children(fine: &VoxelGrid, i: usize, j: usize, k: usize) -> [usize; 8] {
    let mut out = [0; 8];
    for (n, o) in out.iter_mut().enumerate() {
        *o = fine.index(
            2 * i + (n & 1),
            2 * j + ((n >> 1) & 1),
            2 * k + ((n >> 2) & 1),
        );
    }
    out
}

/// The dispersion test of the module doc: whether, for every species of
/// `species`, the counts `count(fine voxel, species)` of every 2x2x2 block
/// of `fine` (one voxel of `coarse`) are consistent with a uniform split at
/// the threshold `z`.
pub(crate) fn blocks_are_uniform(
    fine: &VoxelGrid,
    coarse: &VoxelGrid,
    species: &[usize],
    count: impl Fn(usize, usize) -> i64,
    z: f64,
) -> bool {
    let [nx, ny, nz] = coarse.dims;
    let mut blocks: Vec<[usize; 8]> = Vec::with_capacity(coarse.n_voxels());
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                blocks.push(children(fine, i, j, k));
            }
        }
    }
    for &s in species {
        let (mut excess, mut variance) = (0.0, 0.0);
        // Largest standardized statistic among the populated blocks.
        let (mut worst_block, mut n_tested) = (0.0f64, 0usize);
        for block in &blocks {
            let (mut sum, mut sum_sq) = (0i64, 0i64);
            for &v in block {
                let c = count(v, s);
                sum += c;
                sum_sq += c * c;
            }
            if sum < 1 {
                continue;
            }
            let sf = sum as f64;
            let x = 8.0 * sum_sq as f64 / sf - sf;
            let var = 14.0 * (sf - 1.0) / sf;
            excess += x - 7.0;
            variance += var;
            if sum >= BLOCK_TEST_MIN {
                n_tested += 1;
                worst_block = worst_block.max((x - 7.0) / var.sqrt());
            }
        }
        if variance > 0.0 && excess / variance.sqrt() > z {
            return false;
        }
        if n_tested > 0 && worst_block > z + (2.0 * (n_tested as f64).ln()).sqrt() {
            return false;
        }
    }
    true
}

impl Rdme<'_> {
    /// Turn coarsening on with the dispersion threshold `z` (`None`: off).
    pub fn set_coarsening(&mut self, z: Option<f64>) {
        self.z_max = z;
    }

    /// Whether the counts of every 2x2x2 block are consistent with a
    /// uniform distribution at the threshold `z` (see the module doc).
    /// False if the grid cannot drop a level.
    pub fn is_uniform_in_blocks(&self, z: f64) -> bool {
        if self.level + 1 >= self.setup.levels.len() {
            return false;
        }
        let n_sp = self.setup.n_sp;
        let species: Vec<usize> = (0..n_sp).collect();
        blocks_are_uniform(
            self.grid(),
            &self.setup.levels[self.level + 1].grid,
            &species,
            |v, s| self.counts[v * n_sp + s],
            z,
        )
    }

    /// Merge every 2x2x2 block into one voxel (counts summed). The grid
    /// must not be at its last level.
    pub fn coarsen(&mut self) {
        assert!(
            self.level + 1 < self.setup.levels.len(),
            "rdme: no level to coarsen to"
        );
        let fine = self.grid().clone();
        let coarse = &self.setup.levels[self.level + 1].grid;
        let n_sp = self.setup.n_sp;
        let [nx, ny, nz] = coarse.dims;
        let mut merged = vec![0i64; coarse.n_voxels() * n_sp];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let cv = coarse.index(i, j, k);
                    for v in children(&fine, i, j, k) {
                        for s in 0..n_sp {
                            merged[cv * n_sp + s] += self.counts[v * n_sp + s];
                        }
                    }
                }
            }
        }
        self.counts = merged;
        self.level += 1;
        self.level_history.push((self.t, self.level as u32));
        self.rebuild();
    }

    /// Split every voxel into its 8 children, each molecule going to one of
    /// them at random. The grid must not be at level 0.
    pub fn refine(&mut self) {
        assert!(self.level > 0, "rdme: already on the fine grid");
        let coarse = self.grid().clone();
        let fine = &self.setup.levels[self.level - 1].grid;
        let n_sp = self.setup.n_sp;
        let [nx, ny, nz] = coarse.dims;
        let mut split = vec![0i64; fine.n_voxels() * n_sp];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let cv = coarse.index(i, j, k);
                    let kids = children(fine, i, j, k);
                    for s in 0..n_sp {
                        for _ in 0..self.counts[cv * n_sp + s] {
                            let child = kids[self.rng.below(8) as usize];
                            split[child * n_sp + s] += 1;
                        }
                    }
                }
            }
        }
        self.counts = split;
        self.level -= 1;
        self.level_history.push((self.t, self.level as u32));
        self.rebuild();
    }

    /// Refine down to the fine grid.
    pub fn refine_to_fine(&mut self) {
        while self.level > 0 {
            self.refine();
        }
    }

    /// Drop as many levels as the dispersion test allows. Returns whether
    /// the level changed.
    pub fn coarsen_while_uniform(&mut self) -> bool {
        let Some(z) = self.z_max else {
            return false;
        };
        let before = self.level;
        while self.is_uniform_in_blocks(z) {
            self.coarsen();
        }
        self.level != before
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::config::Config;
    use crate::matrices::SimMatrices;
    use crate::ssa::rdme::RdmeSetup;
    use crate::ssa::rng::SsaRng;

    fn setup(dims: [usize; 3], levels: u32) -> (usize, RdmeSetup) {
        let env = SimEnv::from_config(
            &Config::from_toml_str("[reactions.radiolytic]\ne_aq = 2.6", None).unwrap(),
        )
        .unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new(dims, [10e-9; 3], [0.0; 3]).unwrap();
        let e = env.species.index_of("e_aq").unwrap();
        (
            e,
            RdmeSetup::new(&mat, &grid, vec![5e-9; mat.n_dyn], levels, false).unwrap(),
        )
    }

    #[test]
    fn threshold_scales_with_the_factor() {
        assert_eq!(dispersion_threshold(0.1), 3.0);
        assert_eq!(dispersion_threshold(0.2), 6.0);
    }

    #[test]
    fn a_cluster_blocks_coarsening_and_a_uniform_field_allows_it() {
        let (e, s) = setup([4, 4, 4], 2);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(1, 0));
        // An empty grid is uniform.
        assert!(r.is_uniform_in_blocks(3.0));
        // 40 molecules in one voxel: far from uniform in its block.
        r.inject(&[(0, e, 40)]);
        assert!(!r.is_uniform_in_blocks(3.0));
        // Spread evenly: 5 per voxel of that block.
        let block = children(s.grid(), 0, 0, 0);
        let mut moves = vec![(0, e, -35)];
        moves.extend(block[1..].iter().map(|v| (*v, e, 5)));
        r.inject(&moves);
        assert!(r.is_uniform_in_blocks(3.0));
    }

    #[test]
    fn a_cluster_is_seen_in_a_large_uniform_background() {
        // 512 blocks with a Poisson background of 4 per voxel, and one voxel
        // with 60 more molecules: the aggregate dilutes it, the block test
        // does not.
        let (e, s) = setup([16, 16, 16], 1);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(5, 0));
        let mut y = vec![0.0; s.n_sp];
        y[e] = 4.0 / s.model().molecules_per_umol();
        r.fill_poisson(&y);
        assert!(r.is_uniform_in_blocks(3.0), "a Poisson field is uniform");
        r.inject(&[(1234, e, 60)]);
        assert!(!r.is_uniform_in_blocks(3.0));
    }

    #[test]
    fn coarsening_and_refining_conserve_every_species() {
        let (e, s) = setup([8, 4, 4], 2);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(2, 0));
        let mut y = vec![0.0; s.n_sp];
        y[e] = 3.0 / s.model().molecules_per_umol();
        r.fill_poisson(&y);
        let totals = r.totals();
        let fine_sum: i64 = r.fine_field().iter().sum();
        r.coarsen();
        assert_eq!(r.level(), 1);
        assert_eq!(r.grid().dims, [4, 2, 2]);
        assert_eq!(r.totals(), totals);
        r.coarsen();
        assert_eq!(r.grid().dims, [2, 1, 1]);
        assert_eq!(r.totals(), totals);
        // The fine field spreads each voxel evenly and keeps the total (in
        // units of 1 / 8^2 molecule at every level).
        assert_eq!(r.fine_field().iter().sum::<i64>(), fine_sum);
        assert!(!r.is_uniform_in_blocks(3.0), "last level");
        r.refine_to_fine();
        assert_eq!(r.level(), 0);
        assert_eq!(r.grid().dims, [8, 4, 4]);
        assert_eq!(r.totals(), totals);
        assert_eq!(
            r.level_history.iter().map(|(_, l)| *l).collect::<Vec<_>>(),
            [1, 2, 1, 0]
        );
    }

    #[test]
    fn refining_splits_uniformly() {
        let (e, s) = setup([2, 2, 2], 1);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(3, 0));
        r.inject(&[(0, e, 8000)]);
        r.coarsen();
        assert_eq!(r.voxel(0)[e], 8000);
        r.refine();
        for v in 0..8 {
            // Binomial(8000, 1/8): mean 1000, sd 29.6.
            let n = r.voxel(v)[e];
            assert!((n - 1000).abs() < 5 * 30, "voxel {v}: {n}");
        }
        assert!(r.is_uniform_in_blocks(3.0));
    }

    #[test]
    fn injection_refines_first() {
        let (e, s) = setup([4, 4, 4], 2);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(4, 0));
        r.set_coarsening(Some(3.0));
        assert!(r.coarsen_while_uniform());
        assert_eq!(r.level(), 2);
        r.inject(&[(63, e, 10)]);
        assert_eq!(r.level(), 0);
        assert_eq!(r.voxel(63)[e], 10);
        assert_eq!(r.totals()[e], 10);
    }

    #[test]
    fn the_grid_coarsens_once_diffusion_has_mixed_it() {
        let (e, s) = setup([4, 4, 4], 2);
        let mut r = Rdme::new(&s, 0.0, SsaRng::new(6, 0));
        r.set_coarsening(Some(3.0));
        r.inject(&[(0, e, 3000)]);
        // h^2 / D = 20 ns per voxel: far from mixed after 1 ns.
        r.advance_to(1e-9);
        assert_eq!(r.level(), 0);
        r.advance_to(2e-6);
        assert_eq!(r.level(), 2, "history {:?}", r.level_history);
        assert_eq!(r.totals()[e], 3000);
        let hops_coarse = r.n_hops;
        // On the 1x1x1 grid nothing hops any more.
        r.advance_to(4e-6);
        assert_eq!(r.n_hops, hops_coarse);
    }
}
