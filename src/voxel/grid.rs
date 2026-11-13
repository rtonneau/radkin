//! Uniform Cartesian voxel grid.

use crate::consts::AVOGADRO;
use crate::error::{RadkinError, Result};

/// A uniform Cartesian `nx x ny x nz` grid with per-axis spacing [m].
///
/// Voxel `(i, j, k)` has flat index `i + nx * (j + ny * k)` and covers
/// `[origin + i*dx, origin + (i+1)*dx)` on each axis.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelGrid {
    pub dims: [usize; 3],
    /// Voxel edge lengths [m].
    pub spacing: [f64; 3],
    /// Position of the grid's lower corner [m].
    pub origin: [f64; 3],
}

impl VoxelGrid {
    pub fn new(dims: [usize; 3], spacing: [f64; 3], origin: [f64; 3]) -> Result<Self> {
        if dims.contains(&0) {
            return Err(RadkinError::Config(format!(
                "voxel.grid: every dimension must be >= 1, got {dims:?}"
            )));
        }
        if spacing.iter().any(|s| !(s.is_finite() && *s > 0.0)) {
            return Err(RadkinError::Config(format!(
                "voxel.grid: spacing must be finite and > 0, got {spacing:?}"
            )));
        }
        if origin.iter().any(|o| !o.is_finite()) {
            return Err(RadkinError::Config(format!(
                "voxel.grid: origin must be finite, got {origin:?}"
            )));
        }
        Ok(Self {
            dims,
            spacing,
            origin,
        })
    }

    /// Number of voxels.
    pub fn n_voxels(&self) -> usize {
        self.dims.iter().product()
    }

    /// Volume of one voxel [L].
    pub fn voxel_volume_l(&self) -> f64 {
        self.spacing.iter().product::<f64>() * 1e3
    }

    /// Concentration [mol/L] of one molecule in one voxel.
    pub fn molecule_concentration(&self) -> f64 {
        1.0 / (AVOGADRO * self.voxel_volume_l())
    }

    /// Flat index of voxel `(i, j, k)`.
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.dims[0] * (j + self.dims[1] * k)
    }

    /// `(i, j, k)` of a flat index.
    pub fn coords(&self, idx: usize) -> [usize; 3] {
        let [nx, ny, _] = self.dims;
        [idx % nx, (idx / nx) % ny, idx / (nx * ny)]
    }

    /// Centre of voxel `idx` [m].
    pub fn center(&self, idx: usize) -> [f64; 3] {
        let c = self.coords(idx);
        std::array::from_fn(|a| self.origin[a] + (c[a] as f64 + 0.5) * self.spacing[a])
    }

    /// Flat index of the voxel containing `pos` [m], `None` outside the grid.
    pub fn locate(&self, pos: [f64; 3]) -> Option<usize> {
        let mut c = [0usize; 3];
        for a in 0..3 {
            let f = ((pos[a] - self.origin[a]) / self.spacing[a]).floor();
            if !(f >= 0.0 && f < self.dims[a] as f64) {
                return None;
            }
            c[a] = f as usize;
        }
        Some(self.index(c[0], c[1], c[2]))
    }

    /// Lower and upper corners of the grid [m].
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        let hi = std::array::from_fn(|a| self.origin[a] + self.dims[a] as f64 * self.spacing[a]);
        (self.origin, hi)
    }
}

/// Concentrations [µmol/L, as the kernel] of `n_species` dynamic species on every voxel, stored
/// species-major (structure of arrays): `field(s)[voxel]`.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelState {
    pub n_species: usize,
    pub n_voxels: usize,
    pub data: Vec<f64>,
}

impl VoxelState {
    pub fn zeros(n_species: usize, n_voxels: usize) -> Self {
        Self {
            n_species,
            n_voxels,
            data: vec![0.0; n_species * n_voxels],
        }
    }

    /// Every voxel set to the per-species values `y`.
    pub fn uniform(y: &[f64], n_voxels: usize) -> Self {
        let mut s = Self::zeros(y.len(), n_voxels);
        for (sp, v) in y.iter().enumerate() {
            s.field_mut(sp).fill(*v);
        }
        s
    }

    pub fn field(&self, species: usize) -> &[f64] {
        &self.data[species * self.n_voxels..(species + 1) * self.n_voxels]
    }

    pub fn field_mut(&mut self, species: usize) -> &mut [f64] {
        &mut self.data[species * self.n_voxels..(species + 1) * self.n_voxels]
    }

    /// Copy the species vector of `voxel` into `out`.
    pub fn gather(&self, voxel: usize, out: &mut [f64]) {
        for (sp, o) in out.iter_mut().enumerate() {
            *o = self.data[sp * self.n_voxels + voxel];
        }
    }

    /// Write the species vector `y` into `voxel`.
    pub fn scatter(&mut self, voxel: usize, y: &[f64]) {
        for (sp, v) in y.iter().enumerate() {
            self.data[sp * self.n_voxels + voxel] = *v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> VoxelGrid {
        VoxelGrid::new([3, 4, 5], [1e-9, 2e-9, 3e-9], [-1e-9, 0.0, 0.0]).unwrap()
    }

    #[test]
    fn index_roundtrip() {
        let g = grid();
        assert_eq!(g.n_voxels(), 60);
        for idx in 0..g.n_voxels() {
            let [i, j, k] = g.coords(idx);
            assert_eq!(g.index(i, j, k), idx);
            assert_eq!(g.locate(g.center(idx)), Some(idx));
        }
    }

    #[test]
    fn state_gather_scatter() {
        let mut st = VoxelState::uniform(&[1.0, 2.0], 3);
        st.scatter(1, &[5.0, 6.0]);
        let mut y = [0.0; 2];
        st.gather(1, &mut y);
        assert_eq!(y, [5.0, 6.0]);
        assert_eq!(st.field(0), &[1.0, 5.0, 1.0]);
        assert_eq!(st.field(1), &[2.0, 6.0, 2.0]);
    }

    #[test]
    fn locate_outside() {
        let g = grid();
        assert_eq!(g.locate([-1.5e-9, 0.0, 0.0]), None);
        assert_eq!(g.locate([0.0, 8e-9, 0.0]), None);
        assert_eq!(g.locate([0.0, 0.0, f64::NAN]), None);
        assert_eq!(g.locate([-1e-9, 0.0, 0.0]), Some(0));
    }

    #[test]
    fn rejects_bad_dims_and_spacing() {
        assert!(VoxelGrid::new([0, 1, 1], [1.0; 3], [0.0; 3]).is_err());
        assert!(VoxelGrid::new([1, 1, 1], [1.0, 0.0, 1.0], [0.0; 3]).is_err());
        assert!(VoxelGrid::new([1, 1, 1], [1.0, -1.0, 1.0], [0.0; 3]).is_err());
    }

    #[test]
    fn molecule_concentration_10nm() {
        let g = VoxelGrid::new([1; 3], [10e-9; 3], [0.0; 3]).unwrap();
        // 1 molecule in (10 nm)^3 = 1e-21 L -> 1.66 mM
        assert!((g.molecule_concentration() - 1.0 / (AVOGADRO * 1e-21)).abs() < 1e-15);
    }
}
