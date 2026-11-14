//! Global coarsening primitives of the voxel grid: one level merges every 2x2x2 block.
//!
//! A fine grid with even dims has a coarse grid with half the dims, twice the spacing and the
//! same origin ([`coarse_grid`]). Coarse voxel `(I, J, K)` covers the eight fine voxels
//! `(2I + a, 2J + b, 2K + c)`, `a, b, c` in `{0, 1}`. Concentrations are per volume and the eight
//! children have equal volumes, so the volume average ([`restrict`]) is their plain mean and
//! conserves the amount (concentration * voxel volume) of every species to round-off.
//! Prolongation ([`prolong`]) is piecewise constant, so `restrict(prolong(y)) == y` exactly.
//! [`can_coarsen`] is the per-species spread test that decides whether a level may be dropped.
//!
//! These functions do not depend on the driver and do not allocate, except the `*ed` wrappers,
//! which allocate the output state once per call.

use crate::error::{RadkinError, Result};
use crate::voxel::grid::{VoxelGrid, VoxelState};

/// The grid one coarsening level above `fine`: dims / 2, spacing * 2, same origin.
///
/// Fails if a dimension is odd.
pub fn coarse_grid(fine: &VoxelGrid) -> Result<VoxelGrid> {
    for (axis, n) in ["x", "y", "z"].iter().zip(fine.dims) {
        if n % 2 != 0 {
            return Err(RadkinError::Value(format!(
                "voxel coarsening: grid dims along {axis} ({n}) must be even"
            )));
        }
    }
    VoxelGrid::new(
        fine.dims.map(|n| n / 2),
        fine.spacing.map(|s| 2.0 * s),
        fine.origin,
    )
}

/// Panics unless `grid` has even dims and `state` lives on it.
fn check_fine(grid: &VoxelGrid, state: &VoxelState) {
    assert!(
        grid.dims.iter().all(|n| n % 2 == 0),
        "voxel coarsening: grid dims {:?} must be even",
        grid.dims
    );
    assert_eq!(
        state.n_voxels,
        grid.n_voxels(),
        "voxel coarsening: fine state does not match the grid"
    );
}

/// Panics unless `coarse` matches `fine` one level up.
fn check_pair(grid: &VoxelGrid, fine: &VoxelState, coarse: &VoxelState) {
    check_fine(grid, fine);
    assert_eq!(
        coarse.n_voxels * 8,
        fine.n_voxels,
        "voxel coarsening: coarse state size does not match the fine grid"
    );
    assert_eq!(
        coarse.n_species, fine.n_species,
        "voxel coarsening: species counts differ"
    );
}

/// Mean of eight values, summed pairwise (exact for eight equal values).
fn mean8(v: [f64; 8]) -> f64 {
    0.125 * (((v[0] + v[1]) + (v[2] + v[3])) + ((v[4] + v[5]) + (v[6] + v[7])))
}

/// Flat indices of the eight children of the block whose lowest fine voxel is `base`.
fn children(base: usize, sy: usize, sz: usize) -> [usize; 8] {
    [
        base,
        base + 1,
        base + sy,
        base + sy + 1,
        base + sz,
        base + sz + 1,
        base + sz + sy,
        base + sz + sy + 1,
    ]
}

/// Volume average of every 2x2x2 block of `fine` (on `fine_grid`) into `coarse`.
///
/// # Panics
/// If `fine_grid` has an odd dimension or the state sizes do not match.
pub fn restrict(fine_grid: &VoxelGrid, fine: &VoxelState, coarse: &mut VoxelState) {
    check_pair(fine_grid, fine, coarse);
    let [nx, ny, nz] = fine_grid.dims;
    let [cx, cy, cz] = [nx / 2, ny / 2, nz / 2];
    let (sy, sz) = (nx, nx * ny);
    for s in 0..fine.n_species {
        let f = fine.field(s);
        let c = coarse.field_mut(s);
        for k in 0..cz {
            for j in 0..cy {
                let row = cx * (j + cy * k);
                let base = nx * (2 * j + ny * 2 * k);
                for i in 0..cx {
                    c[row + i] = mean8(children(base + 2 * i, sy, sz).map(|idx| f[idx]));
                }
            }
        }
    }
}

/// [`restrict`] into a newly allocated coarse state.
pub fn restricted(fine_grid: &VoxelGrid, fine: &VoxelState) -> VoxelState {
    let mut coarse = VoxelState::zeros(fine.n_species, fine.n_voxels / 8);
    restrict(fine_grid, fine, &mut coarse);
    coarse
}

/// Piecewise-constant prolongation: every fine voxel of `fine_grid` takes the value of its
/// coarse parent in `coarse`.
///
/// # Panics
/// If `fine_grid` has an odd dimension or the state sizes do not match.
pub fn prolong(fine_grid: &VoxelGrid, coarse: &VoxelState, fine: &mut VoxelState) {
    check_pair(fine_grid, fine, coarse);
    let [nx, ny, nz] = fine_grid.dims;
    let [cx, cy, _] = [nx / 2, ny / 2, nz / 2];
    for s in 0..fine.n_species {
        let c = coarse.field(s);
        let f = fine.field_mut(s);
        for k in 0..nz {
            for j in 0..ny {
                let crow = &c[cx * (j / 2 + cy * (k / 2))..][..cx];
                let frow = &mut f[nx * (j + ny * k)..][..nx];
                for (i, v) in frow.iter_mut().enumerate() {
                    *v = crow[i / 2];
                }
            }
        }
    }
}

/// [`prolong`] into a newly allocated fine state.
pub fn prolonged(fine_grid: &VoxelGrid, coarse: &VoxelState) -> VoxelState {
    let mut fine = VoxelState::zeros(coarse.n_species, fine_grid.n_voxels());
    prolong(fine_grid, coarse, &mut fine);
    fine
}

/// Widen `[lo, hi]` with the four face-neighbour values at `start + {0, u, v, u + v}`.
fn extend4(f: &[f64], start: usize, u: usize, v: usize, lo: &mut f64, hi: &mut f64) {
    for idx in [start, start + u, start + v, start + u + v] {
        *lo = lo.min(f[idx]);
        *hi = hi.max(f[idx]);
    }
}

/// Whether `fine` (on `fine_grid`) may drop one coarsening level.
///
/// For every species and every 2x2x2 block, the spread (max - min) over the eight children and
/// their face neighbours must not exceed `factor * (atol + rtol * |mean|)`, with `mean` the mean
/// of the eight children. Face neighbours outside the domain are skipped (zero-flux walls). A
/// non-finite value fails the test. `atol`, `rtol` and `factor` are expected to be >= 0.
///
/// # Panics
/// If `fine_grid` has an odd dimension or `fine` does not live on it.
pub fn can_coarsen(
    fine_grid: &VoxelGrid,
    fine: &VoxelState,
    atol: f64,
    rtol: f64,
    factor: f64,
) -> bool {
    check_fine(fine_grid, fine);
    let [nx, ny, nz] = fine_grid.dims;
    let (sy, sz) = (nx, nx * ny);
    for s in 0..fine.n_species {
        let f = fine.field(s);
        for k in (0..nz).step_by(2) {
            for j in (0..ny).step_by(2) {
                for i in (0..nx).step_by(2) {
                    let base = i + nx * (j + ny * k);
                    let vals = children(base, sy, sz).map(|idx| f[idx]);
                    let mean = mean8(vals);
                    let mut lo = f64::INFINITY;
                    let mut hi = f64::NEG_INFINITY;
                    for v in vals {
                        lo = lo.min(v);
                        hi = hi.max(v);
                    }
                    if i > 0 {
                        extend4(f, base - 1, sy, sz, &mut lo, &mut hi);
                    }
                    if i + 2 < nx {
                        extend4(f, base + 2, sy, sz, &mut lo, &mut hi);
                    }
                    if j > 0 {
                        extend4(f, base - sy, 1, sz, &mut lo, &mut hi);
                    }
                    if j + 2 < ny {
                        extend4(f, base + 2 * sy, 1, sz, &mut lo, &mut hi);
                    }
                    if k > 0 {
                        extend4(f, base - sz, 1, sy, &mut lo, &mut hi);
                    }
                    if k + 2 < nz {
                        extend4(f, base + 2 * sz, 1, sy, &mut lo, &mut hi);
                    }
                    // Every voxel is a child of some block, so a NaN anywhere makes some
                    // block's mean non-finite.
                    let tol = factor * (atol + rtol * mean.abs());
                    if !(mean.is_finite() && hi - lo <= tol) {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    const ATOL: f64 = 1e-6;
    const RTOL: f64 = 1e-3;
    const FACTOR: f64 = 0.1;

    fn fine_grid() -> VoxelGrid {
        VoxelGrid::new([4, 6, 2], [1e-9, 2e-9, 3e-9], [-1e-9, 5e-9, 0.0]).unwrap()
    }

    /// Deterministic pseudo-random values in [0.5, 10.5) with mixed magnitudes.
    fn lcg_state(n_species: usize, n_voxels: usize, seed: u64) -> VoxelState {
        let mut x = seed;
        let mut st = VoxelState::zeros(n_species, n_voxels);
        for v in st.data.iter_mut() {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *v = 0.5 + 10.0 * ((x >> 11) as f64 / (1u64 << 53) as f64);
        }
        st
    }

    fn totals(grid: &VoxelGrid, st: &VoxelState) -> Vec<f64> {
        (0..st.n_species)
            .map(|s| st.field(s).iter().sum::<f64>() * grid.voxel_volume_l())
            .collect()
    }

    #[test]
    fn coarse_grid_halves_dims_doubles_spacing() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        assert_eq!(c.dims, [2, 3, 1]);
        assert_eq!(c.spacing, [2e-9, 4e-9, 6e-9]);
        assert_eq!(c.origin, g.origin);
        assert_eq!(c.bounds(), g.bounds());
        assert_eq!(c.voxel_volume_l(), 8.0 * g.voxel_volume_l());
        let odd = VoxelGrid::new([4, 3, 2], [1e-9; 3], [0.0; 3]).unwrap();
        let err = coarse_grid(&odd).unwrap_err().to_string();
        assert!(err.contains("along y (3)"), "{err}");
        assert!(coarse_grid(&c).is_err());
    }

    #[test]
    fn restrict_conserves_totals() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let fine = lcg_state(3, g.n_voxels(), 7);
        let coarse = restricted(&g, &fine);
        assert_eq!(coarse.n_voxels, c.n_voxels());
        for (tf, tc) in totals(&g, &fine).iter().zip(totals(&c, &coarse)) {
            assert!((tf - tc).abs() <= 1e-14 * tf.abs(), "{tf} vs {tc}");
        }
    }

    #[test]
    fn restrict_is_block_mean() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let fine = lcg_state(2, g.n_voxels(), 3);
        let coarse = restricted(&g, &fine);
        for s in 0..2 {
            for idx in 0..c.n_voxels() {
                let [ci, cj, ck] = c.coords(idx);
                let mut sum = 0.0;
                for (a, b, d) in (0..8).map(|n| (n & 1, (n >> 1) & 1, n >> 2)) {
                    sum += fine.field(s)[g.index(2 * ci + a, 2 * cj + b, 2 * ck + d)];
                }
                let want = sum / 8.0;
                let got = coarse.field(s)[idx];
                assert!((got - want).abs() <= 1e-15 * want.abs(), "{got} vs {want}");
            }
        }
    }

    #[test]
    fn prolong_of_restrict_is_exact_for_blockwise_constant() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let y = lcg_state(3, c.n_voxels(), 11);
        let mut x = VoxelState::zeros(3, g.n_voxels());
        for s in 0..3 {
            for idx in 0..g.n_voxels() {
                let [i, j, k] = g.coords(idx);
                x.field_mut(s)[idx] = y.field(s)[c.index(i / 2, j / 2, k / 2)];
            }
        }
        let r = restricted(&g, &x);
        assert_eq!(r, y);
        assert_eq!(prolonged(&g, &r), x);
    }

    #[test]
    fn restrict_of_prolong_is_identity() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let y = lcg_state(2, c.n_voxels(), 5);
        assert_eq!(restricted(&g, &prolonged(&g, &y)), y);
    }

    #[test]
    fn in_place_buffers_are_overwritten() {
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let fine = lcg_state(2, g.n_voxels(), 13);
        let mut coarse = VoxelState::uniform(&[f64::NAN, -1.0], c.n_voxels());
        restrict(&g, &fine, &mut coarse);
        assert_eq!(coarse, restricted(&g, &fine));
        let mut back = VoxelState::uniform(&[f64::NAN, -1.0], g.n_voxels());
        prolong(&g, &coarse, &mut back);
        assert_eq!(back, prolonged(&g, &coarse));
    }

    #[test]
    fn uniform_field_can_coarsen() {
        let g = fine_grid();
        let st = VoxelState::uniform(&[0.0, 3.5, 1e3], g.n_voxels());
        assert!(can_coarsen(&g, &st, ATOL, RTOL, FACTOR));
        // a zero field coarsens even with zero tolerances
        let zero = VoxelState::zeros(1, g.n_voxels());
        assert!(can_coarsen(&g, &zero, 0.0, 0.0, FACTOR));
    }

    #[test]
    fn small_noise_below_tolerance_can_coarsen() {
        let g = fine_grid();
        let mut st = VoxelState::uniform(&[100.0], g.n_voxels());
        // tol = 0.1 * (1e-6 + 1e-3 * 100) ~ 1e-2
        for (n, v) in st.field_mut(0).iter_mut().enumerate() {
            *v += if n % 3 == 0 { 4e-3 } else { -4e-3 };
        }
        assert!(can_coarsen(&g, &st, ATOL, RTOL, FACTOR));
        st.field_mut(0)[5] += 2e-2;
        assert!(!can_coarsen(&g, &st, ATOL, RTOL, FACTOR));
    }

    #[test]
    fn single_spike_cannot_coarsen() {
        let g = fine_grid();
        for idx in [0, 9, g.n_voxels() - 1] {
            let mut st = VoxelState::uniform(&[1.0, 2.0], g.n_voxels());
            st.field_mut(1)[idx] = 10.0;
            assert!(!can_coarsen(&g, &st, ATOL, RTOL, FACTOR), "spike at {idx}");
        }
    }

    #[test]
    fn face_neighbour_blocks_coarsening() {
        // Each block is internally uniform, but adjacent blocks differ: the spread over the
        // face neighbours must catch it.
        let g = fine_grid();
        let c = coarse_grid(&g).unwrap();
        let mut y = VoxelState::uniform(&[1.0], c.n_voxels());
        y.field_mut(0)[c.index(1, 2, 0)] = 1.5;
        let x = prolonged(&g, &y);
        assert!(!can_coarsen(&g, &x, ATOL, RTOL, FACTOR));
        // ... unless the jump is within tolerance
        y.field_mut(0)[c.index(1, 2, 0)] = 1.0 + 5e-5;
        assert!(can_coarsen(&g, &prolonged(&g, &y), ATOL, RTOL, FACTOR));
    }

    #[test]
    fn single_block_grid_skips_walls() {
        let g = VoxelGrid::new([2, 2, 2], [1e-9; 3], [0.0; 3]).unwrap();
        assert!(can_coarsen(
            &g,
            &VoxelState::uniform(&[2.0], 8),
            ATOL,
            RTOL,
            FACTOR
        ));
        let mut st = VoxelState::uniform(&[2.0], 8);
        st.field_mut(0)[7] = 3.0;
        assert!(!can_coarsen(&g, &st, ATOL, RTOL, FACTOR));
    }

    #[test]
    fn non_finite_cannot_coarsen() {
        let g = fine_grid();
        for bad in [f64::NAN, f64::INFINITY] {
            let mut st = VoxelState::uniform(&[1.0], g.n_voxels());
            st.field_mut(0)[17] = bad;
            assert!(!can_coarsen(&g, &st, ATOL, RTOL, FACTOR), "{bad}");
        }
    }
}
