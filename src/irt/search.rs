//! Pair search of the IRT stage: every pair of particles closer than a
//! cutoff radius, found with a uniform grid of cells.

use crate::irt::sample::Vec3;

/// No particle (end of a cell list).
const NONE: u32 = u32::MAX;
/// Below this many particles every pair is tested directly.
const BRUTE_FORCE: usize = 64;
/// Most cells per axis.
const MAX_CELLS: usize = 256;

/// Scratch of the pair search, reused between calls.
#[derive(Debug, Clone, Default)]
pub struct CellGrid {
    /// First particle of each cell.
    heads: Vec<u32>,
    /// Next particle of the same cell.
    next: Vec<u32>,
    /// Cell coordinates of each particle.
    cell: Vec<[u32; 3]>,
}

fn dist2(a: Vec3, b: Vec3) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

impl CellGrid {
    pub fn new() -> Self {
        Self::default()
    }

    /// Call `f(i, j, distance^2)` with `i < j` for every pair of `pos` closer
    /// than `cutoff`, in an order that depends on the positions only.
    pub fn for_each_pair(
        &mut self,
        pos: &[Vec3],
        cutoff: f64,
        mut f: impl FnMut(usize, usize, f64),
    ) {
        let n = pos.len();
        let c2 = cutoff * cutoff;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in pos {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        // Cells at least `cutoff` wide, so a pair is in the same cell or in
        // two neighbouring ones.
        let dims: [usize; 3] = std::array::from_fn(|a| {
            let cells = ((hi[a] - lo[a]) / cutoff).floor();
            if cells.is_finite() && cells >= 1.0 {
                (cells as usize).min(MAX_CELLS)
            } else {
                1
            }
        });
        if n <= BRUTE_FORCE || dims.iter().all(|d| *d <= 2) {
            for i in 0..n {
                for j in i + 1..n {
                    let d2 = dist2(pos[i], pos[j]);
                    if d2 < c2 {
                        f(i, j, d2);
                    }
                }
            }
            return;
        }
        let width: [f64; 3] = std::array::from_fn(|a| (hi[a] - lo[a]) / dims[a] as f64);
        let n_cells = dims[0] * dims[1] * dims[2];
        self.heads.clear();
        self.heads.resize(n_cells, NONE);
        self.next.clear();
        self.next.resize(n, NONE);
        self.cell.clear();
        let flat =
            |c: [u32; 3]| c[0] as usize + dims[0] * (c[1] as usize + dims[1] * c[2] as usize);
        for p in pos {
            let c: [u32; 3] = std::array::from_fn(|a| {
                (((p[a] - lo[a]) / width[a]) as usize).min(dims[a] - 1) as u32
            });
            self.cell.push(c);
        }
        // Filled from the last particle, so each list is in rising order.
        for i in (0..n).rev() {
            let k = flat(self.cell[i]);
            self.next[i] = self.heads[k];
            self.heads[k] = i as u32;
        }
        for i in 0..n {
            let ci = self.cell[i];
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let c = [ci[0] as i64 + dx, ci[1] as i64 + dy, ci[2] as i64 + dz];
                        if (0..3).any(|a| c[a] < 0 || c[a] >= dims[a] as i64) {
                            continue;
                        }
                        let mut j = self.heads[flat([c[0] as u32, c[1] as u32, c[2] as u32])];
                        while j != NONE {
                            let ju = j as usize;
                            if ju > i {
                                let d2 = dist2(pos[i], pos[ju]);
                                if d2 < c2 {
                                    f(i, ju, d2);
                                }
                            }
                            j = self.next[ju];
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssa::rng::SsaRng;

    fn pairs(pos: &[Vec3], cutoff: f64) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        CellGrid::new().for_each_pair(pos, cutoff, |i, j, d2| {
            assert!(i < j);
            assert!((d2 - dist2(pos[i], pos[j])).abs() < 1e-30);
            out.push((i, j));
        });
        out.sort_unstable();
        out
    }

    fn brute(pos: &[Vec3], cutoff: f64) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for i in 0..pos.len() {
            for j in i + 1..pos.len() {
                if dist2(pos[i], pos[j]) < cutoff * cutoff {
                    out.push((i, j));
                }
            }
        }
        out
    }

    #[test]
    fn grid_search_finds_the_pairs_of_the_direct_search() {
        let mut rng = SsaRng::new(3, 0);
        // an elongated cloud, like a track: 1 x 1 x 20 um
        let pos: Vec<Vec3> = (0..3000)
            .map(|_| {
                [
                    1e-6 * rng.uniform(),
                    1e-6 * rng.uniform(),
                    20e-6 * rng.uniform(),
                ]
            })
            .collect();
        for cutoff in [30e-9, 200e-9, 0.9e-6, 5e-6, 50e-6] {
            let got = pairs(&pos, cutoff);
            let want = brute(&pos, cutoff);
            assert!(!want.is_empty());
            assert_eq!(got, want, "cutoff {cutoff}");
            // no pair twice
            let mut dedup = got.clone();
            dedup.dedup();
            assert_eq!(dedup.len(), got.len());
        }
    }

    #[test]
    fn small_and_degenerate_sets() {
        assert!(pairs(&[], 1.0).is_empty());
        assert!(pairs(&[[0.0; 3]], 1.0).is_empty());
        // identical positions
        let pos = vec![[1.0, 2.0, 3.0]; 100];
        assert_eq!(pairs(&pos, 1e-9).len(), 100 * 99 / 2);
        // a plane: no extent along z
        let mut rng = SsaRng::new(4, 0);
        let pos: Vec<Vec3> = (0..500)
            .map(|_| [rng.uniform(), rng.uniform(), 0.5])
            .collect();
        assert_eq!(pairs(&pos, 0.05), brute(&pos, 0.05));
    }
}
