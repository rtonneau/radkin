//! Implicit diffusion step on the voxel grid with zero-flux boundaries.
//!
//! The 3D operator `D * Laplacian` is split into the 1D second differences
//! `Ax`, `Ay`, `Az` (ghost-cell mirrored, zero-flux). On a uniform grid with
//! Neumann boundaries they commute, so the product of 1D theta-scheme factors
//! `(I - theta*dt*A)^-1 (I + (1-theta)*dt*A)` over the three axes is an ADI
//! (locally one-dimensional) scheme: second order for Crank-Nicolson
//! (`theta = 1/2`), unconditionally stable, and mass-conserving to round-off.
//! Backward Euler (`theta = 1`) factors are inverses of M-matrices, so they keep
//! concentrations non-negative; they damp the Crank-Nicolson oscillations of
//! spiky data (Rannacher start, [`Diffusion::step_rannacher`]).

use crate::voxel::grid::{VoxelGrid, VoxelState};

/// Time discretisation of one diffusion step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheme {
    CrankNicolson,
    BackwardEuler,
}

impl Scheme {
    fn theta(self) -> f64 {
        match self {
            Self::CrankNicolson => 0.5,
            Self::BackwardEuler => 1.0,
        }
    }
}

/// Diffusion operator of a grid, with per-species coefficients [m^2/s].
#[derive(Debug, Clone)]
pub struct Diffusion {
    grid: VoxelGrid,
    coeffs: Vec<f64>,
    line: Vec<f64>,
    rhs: Vec<f64>,
    cprime: Vec<f64>,
}

impl Diffusion {
    /// `coeffs[s]` is the diffusion coefficient of species `s` (0 = immobile).
    pub fn new(grid: VoxelGrid, coeffs: Vec<f64>) -> Self {
        let n = grid.dims.iter().copied().max().unwrap_or(1);
        Self {
            grid,
            coeffs,
            line: vec![0.0; n],
            rhs: vec![0.0; n],
            cprime: vec![0.0; n],
        }
    }

    pub fn grid(&self) -> &VoxelGrid {
        &self.grid
    }

    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// Advance every mobile species of `state` by `dt` [s].
    pub fn step(&mut self, state: &mut VoxelState, dt: f64, scheme: Scheme) {
        debug_assert_eq!(state.n_species, self.coeffs.len());
        debug_assert_eq!(state.n_voxels, self.grid.n_voxels());
        for sp in 0..self.coeffs.len() {
            let d = self.coeffs[sp];
            if d == 0.0 || dt == 0.0 {
                continue;
            }
            for axis in 0..3 {
                if self.grid.dims[axis] < 2 {
                    continue;
                }
                let r = d * dt / (self.grid.spacing[axis] * self.grid.spacing[axis]);
                self.axis_pass(state.field_mut(sp), axis, r, scheme.theta());
            }
        }
    }

    /// `n_be` backward Euler steps of `dt / n_be` (Rannacher start).
    pub fn step_rannacher(&mut self, state: &mut VoxelState, dt: f64, n_be: usize) {
        let n = n_be.max(1);
        for _ in 0..n {
            self.step(state, dt / n as f64, Scheme::BackwardEuler);
        }
    }

    fn axis_pass(&mut self, u: &mut [f64], axis: usize, r: f64, theta: f64) {
        let [nx, ny, nz] = self.grid.dims;
        // (line length, stride along the line, outer counts and strides).
        let (n, stride, outer_a, outer_b, sa, sb) = match axis {
            0 => (nx, 1, ny, nz, nx, nx * ny),
            1 => (ny, nx, nx, nz, 1, nx * ny),
            _ => (nz, nx * ny, nx, ny, 1, nx),
        };
        let ri = theta * r;
        let re = (1.0 - theta) * r;
        for b in 0..outer_b {
            for a in 0..outer_a {
                let start = a * sa + b * sb;
                for i in 0..n {
                    self.line[i] = u[start + i * stride];
                }
                explicit_part(&self.line[..n], re, &mut self.rhs[..n]);
                thomas_neumann(&mut self.rhs[..n], ri, &mut self.cprime[..n]);
                for i in 0..n {
                    u[start + i * stride] = self.rhs[i];
                }
            }
        }
    }
}

/// `out = (I + re * A) u` with `A` the zero-flux second difference (n >= 2).
fn explicit_part(u: &[f64], re: f64, out: &mut [f64]) {
    let n = u.len();
    if re == 0.0 {
        out.copy_from_slice(u);
        return;
    }
    out[0] = u[0] + re * (u[1] - u[0]);
    for i in 1..n - 1 {
        out[i] = u[i] + re * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
    }
    out[n - 1] = u[n - 1] + re * (u[n - 2] - u[n - 1]);
}

/// Solve `(I - ri * A) x = d` in place (Thomas), `A` zero-flux, n >= 2.
fn thomas_neumann(d: &mut [f64], ri: f64, cp: &mut [f64]) {
    let n = d.len();
    let off = -ri;
    let diag = |i: usize| {
        if i == 0 || i == n - 1 {
            1.0 + ri
        } else {
            1.0 + 2.0 * ri
        }
    };
    let mut denom = diag(0);
    cp[0] = off / denom;
    d[0] /= denom;
    for i in 1..n {
        denom = diag(i) - off * cp[i - 1];
        cp[i] = off / denom;
        d[i] = (d[i] - off * d[i - 1]) / denom;
    }
    for i in (0..n - 1).rev() {
        d[i] -= cp[i] * d[i + 1];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const NM: f64 = 1e-9;

    fn gaussian(grid: &VoxelGrid, sigma2: f64, dim: usize) -> Vec<f64> {
        let (lo, hi) = grid.bounds();
        let c: [f64; 3] = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
        (0..grid.n_voxels())
            .map(|i| {
                let p = grid.center(i);
                let r2: f64 = (0..3).map(|a| (p[a] - c[a]).powi(2)).sum();
                (2.0 * PI * sigma2).powf(-(dim as f64) / 2.0) * (-r2 / (2.0 * sigma2)).exp()
            })
            .collect()
    }

    fn run_gaussian(dims: [usize; 3], dim: usize) -> f64 {
        let grid = VoxelGrid::new(dims, [NM; 3], [0.0; 3]).unwrap();
        let d = 1e-9;
        let (s0, t) = ((4.0 * NM).powi(2), 4.0e-9);
        let mut st = VoxelState::zeros(1, grid.n_voxels());
        st.field_mut(0).copy_from_slice(&gaussian(&grid, s0, dim));
        let mut diff = Diffusion::new(grid.clone(), vec![d]);
        let n = 200;
        for _ in 0..n {
            diff.step(&mut st, t / n as f64, Scheme::CrankNicolson);
        }
        let exact = gaussian(&grid, s0 + 2.0 * d * t, dim);
        let peak = exact.iter().cloned().fold(0.0, f64::max);
        st.field(0)
            .iter()
            .zip(&exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / peak
    }

    #[test]
    fn gaussian_3d() {
        let err = run_gaussian([61, 61, 61], 3);
        assert!(err < 1e-2, "max rel error {err}");
    }

    #[test]
    fn gaussian_1d_and_2d_degenerate() {
        assert!(run_gaussian([81, 1, 1], 1) < 1e-2);
        assert!(run_gaussian([61, 61, 1], 2) < 1e-2);
    }

    #[test]
    fn mass_conserved_and_be_positive() {
        let grid = VoxelGrid::new([15, 9, 7], [2.0 * NM, NM, 3.0 * NM], [0.0; 3]).unwrap();
        let mut st = VoxelState::zeros(2, grid.n_voxels());
        st.field_mut(0)[grid.index(0, 4, 6)] = 1.0;
        st.field_mut(1)[grid.index(7, 0, 3)] = 1.0;
        let mut diff = Diffusion::new(grid, vec![5e-9, 0.0]);
        // Large step: D dt / dx^2 = 50 on the finest axis.
        diff.step_rannacher(&mut st, 1e-8, 4);
        assert!(st.field(0).iter().all(|v| *v >= 0.0));
        let m0: f64 = st.field(0).iter().sum();
        assert!((m0 - 1.0).abs() < 1e-12, "{m0}");
        for _ in 0..20 {
            diff.step(&mut st, 1e-8, Scheme::CrankNicolson);
        }
        let m0: f64 = st.field(0).iter().sum();
        assert!((m0 - 1.0).abs() < 1e-12, "{m0}");
        // Immobile species untouched.
        assert_eq!(st.field(1).iter().sum::<f64>(), 1.0);
        assert_eq!(st.field(1).iter().filter(|v| **v != 0.0).count(), 1);
    }

    #[test]
    fn second_order_in_time() {
        let grid = VoxelGrid::new([41, 1, 1], [NM; 3], [0.0; 3]).unwrap();
        let init = gaussian(&grid, (3.0 * NM).powi(2), 1);
        let run = |n: usize| {
            let mut st = VoxelState::zeros(1, grid.n_voxels());
            st.field_mut(0).copy_from_slice(&init);
            let mut diff = Diffusion::new(grid.clone(), vec![1e-9]);
            for _ in 0..n {
                diff.step(&mut st, 4e-9 / n as f64, Scheme::CrankNicolson);
            }
            st.field(0).to_vec()
        };
        let reference = run(1024);
        let e = |n: usize| {
            run(n)
                .iter()
                .zip(&reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        };
        let ratio = e(8) / e(16);
        assert!((3.5..4.5).contains(&ratio), "ratio {ratio}");
    }
}
