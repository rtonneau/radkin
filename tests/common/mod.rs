//! The seedless pH 7 background of `tests/ssa_rdme.rs` and
//! `tests/ssa_rdme_multigrid.rs` (ADR 0006): the RDME with the mesh
//! correction against the unit-cell SSA of the same box.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::matrices::SimMatrices;
use radkin::ssa::{Cell, RdmeSetup, SsaModel, SsaRng};
use radkin::voxel::config::DEFAULT_DIFFUSION;
use radkin::voxel::grid::VoxelGrid;

/// Allowed distance in standard errors.
pub const K: f64 = 4.0;

/// Run time of the background [s].
pub const T_END: f64 = 1e-5;

/// Voxels per box edge.
const N: usize = 32;

/// Trajectories of the unit-cell reference.
const CELL_RUNS: usize = 20_000;

/// `f(i)` for `i` in `0..n` on every core, results in index order.
pub fn par_map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let threads = std::thread::available_parallelism().map_or(1, |t| t.get());
    let chunk = n.div_ceil(threads).max(1);
    std::thread::scope(|scope| {
        let f = &f;
        let handles: Vec<_> = (0..n)
            .step_by(chunk)
            .map(|lo| scope.spawn(move || (lo..(lo + chunk).min(n)).map(f).collect::<Vec<T>>()))
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    })
}

/// Mean and standard error of every column of `rows`.
fn mean_sem(rows: &[[f64; 2]]) -> [(f64, f64); 2] {
    let n = rows.len() as f64;
    [0, 1].map(|c| {
        let mean = rows.iter().map(|r| r[c]).sum::<f64>() / n;
        let var = rows.iter().map(|r| (r[c] - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (mean, (var / n).sqrt())
    })
}

/// CONFIG_PAPER_01 at pH 7 in a seedless box of `32^3` voxels, with the
/// default diffusion coefficients, except O2: nothing reacts with it here,
/// and its hops would cost more than the ions'. Only water autoionisation
/// `H2O <-> H3O_plus + OH_minus` runs: its backward direction is
/// bimolecular (h* = 3.4 nm), its forward one first order.
pub struct Background {
    pub env: SimEnv,
    /// Mesh-corrected setup with `depth` coarsening levels.
    pub setup: RdmeSetup,
    /// H3O_plus and OH_minus.
    pub ions: [usize; 2],
    /// Initial state [µmol/L].
    pub y0: Vec<f64>,
}

impl Background {
    pub fn new(h_nm: f64, depth: u32) -> Self {
        let env = SimEnv::from_config(&Config::load("data/CONFIG_PAPER_01.toml").unwrap()).unwrap();
        let mat = SimMatrices::build(&env).unwrap();
        let grid = VoxelGrid::new([N; 3], [h_nm * 1e-9; 3], [0.0; 3]).unwrap();
        let diffusion = env
            .species
            .sorted_dyn_labels()
            .iter()
            .map(|s| match DEFAULT_DIFFUSION.iter().find(|(x, _)| x == s) {
                Some((x, d)) if *x != "O2" => *d,
                _ => 0.0,
            })
            .collect();
        let setup = RdmeSetup::new(&mat, &grid, diffusion, depth, true).unwrap();
        let ions = ["H3O_plus", "OH_minus"].map(|s| env.species.index_of(s).unwrap());
        let y0 = env.initial_values();
        Self {
            env,
            setup,
            ions,
            y0,
        }
    }

    /// Change of the ions between two count vectors.
    pub fn change(&self, start: &[i64], end: &[i64]) -> [f64; 2] {
        self.ions.map(|s| (end[s] - start[s]) as f64)
    }

    /// Changes of the ions over [`T_END`] in the unit-cell SSA of the box
    /// volume, from Poisson starts.
    fn unit_cell(&self) -> Vec<[f64; 2]> {
        let mat = SimMatrices::build(&self.env).unwrap();
        let grid = self.setup.grid();
        let model = SsaModel::new(&mat, grid.voxel_volume_l() * grid.n_voxels() as f64).unwrap();
        let k = model.molecules_per_umol();
        par_map(CELL_RUNS, |i| {
            let mut rng = SsaRng::new(77, i as u64);
            let start: Vec<i64> = self.y0.iter().map(|y| rng.poisson(y * k)).collect();
            let mut cell = Cell::new(&model, start.clone(), 0.0);
            cell.advance(T_END, 0.0, &mut rng);
            self.change(&start, &cell.counts)
        })
    }

    /// Distance in standard errors of the mean RDME change of each ion from
    /// the unit-cell one; each must be within [`K`].
    pub fn compare(&self, name: &str, rdme: &[[f64; 2]]) -> [f64; 2] {
        let h_nm = self.setup.grid().spacing[0] * 1e9;
        let start = self.y0[self.ions[0]] * self.setup.model().molecules_per_umol();
        let (r, c) = (mean_sem(rdme), mean_sem(&self.unit_cell()));
        [0, 1].map(|i| {
            let ((mr, sr), (mc, sc)) = (r[i], c[i]);
            let z = (mr - mc) / (sr * sr + sc * sc).sqrt();
            let line = format!(
                "{name}, {h_nm} nm, ion {i}: {mr:.4} +- {sr:.4} (rdme, {} trajectories) vs \
                 {mc:.4} +- {sc:.4} (unit cell) of {:.2} molecules per box: {z:.2} SEM",
                rdme.len(),
                start * self.setup.grid().n_voxels() as f64
            );
            println!("{line}");
            assert!(z.abs() <= K, "{line}");
            z
        })
    }
}

/// The pooled distance of the sizes, `sum z / sqrt(n)`, is within [`K`]:
/// a drift too small to fail one size fails the pool.
pub fn assert_pooled(name: &str, z: &[[f64; 2]]) {
    for i in 0..2 {
        let pooled = z.iter().map(|x| x[i]).sum::<f64>() / (z.len() as f64).sqrt();
        println!("{name}, ion {i}: pooled {pooled:.2} SEM");
        assert!(pooled.abs() <= K, "{name}, ion {i}: pooled {pooled:.2} SEM");
    }
}
