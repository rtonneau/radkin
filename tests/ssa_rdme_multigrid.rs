//! Multi-grid RDME (`[ssa.species_levels]`, ADR 0005): reactions between
//! species on different levels of the voxel hierarchy.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::consts::AVOGADRO;
use radkin::matrices::SimMatrices;
use radkin::ssa::{Multigrid, MultigridSetup, RdmeSetup, SsaRng};
use radkin::voxel::grid::VoxelGrid;

/// Environment, multi-grid setup and species index lookup for `toml` on a
/// cubic grid of spacing `h`, with per-species diffusion and levels.
fn setup(
    toml: &str,
    dims: [usize; 3],
    h: f64,
    depth: u32,
    correct: bool,
    species: &[(&str, f64, u32)],
) -> (SimEnv, MultigridSetup) {
    let env = SimEnv::from_config(&Config::from_toml_str(toml, None).unwrap()).unwrap();
    let mat = SimMatrices::build(&env).unwrap();
    let grid = VoxelGrid::new(dims, [h; 3], [0.0; 3]).unwrap();
    let (mut d, mut l) = (vec![0.0; mat.n_dyn], vec![0; mat.n_dyn]);
    for (name, coeff, level) in species {
        let s = env.species.index_of(name).unwrap();
        d[s] = *coeff;
        l[s] = *level;
    }
    let base = RdmeSetup::new(&mat, &grid, d, depth, correct).unwrap();
    (env, MultigridSetup::new(base, &l).unwrap())
}

/// `out += rate * (discrete Laplacian of m)` on `dims` with zero-flux faces.
fn add_laplacian(dims: [usize; 3], rate: f64, m: &[f64], out: &mut [f64]) {
    let idx = |c: [usize; 3]| c[0] + dims[0] * (c[1] + dims[1] * c[2]);
    for k in 0..dims[2] {
        for j in 0..dims[1] {
            for i in 0..dims[0] {
                let c = [i, j, k];
                let v = idx(c);
                for a in 0..3 {
                    for up in [false, true] {
                        let mut n = c;
                        if up {
                            if c[a] + 1 == dims[a] {
                                continue;
                            }
                            n[a] += 1;
                        } else {
                            if c[a] == 0 {
                                continue;
                            }
                            n[a] -= 1;
                        }
                        out[v] += rate * (m[idx(n)] - m[v]);
                    }
                }
            }
        }
    }
}

#[test]
fn ssa_multigrid_conversion_to_a_finer_level_matches_the_mean_field() {
    // A (H_r) on level 1, converted at rate k into B (H2) on level 0, which
    // lands uniformly in one of the 8 fine voxels of A's voxel (Erban and
    // Winkelmann 2024, eq. 20). The means of the mgRDME obey (eq. 21-22)
    //   dA/dt = D_A/h_A^2 L_A A - k A,
    //   dB/dt = D_B/h_B^2 L_B B + k/8 sum_{children} A,
    // integrated here with RK4.
    let toml = "[[reactions.k_reaction]]\nreaction = \"H_r -> H2\"\nk_value = 2.0e6";
    let (h, d_a, d_b, k) = (50e-9, 5e-9, 2e-9, 2.0e6);
    let (env, s) = setup(
        toml,
        [8, 2, 2],
        h,
        1,
        false,
        &[("H_r", d_a, 1), ("H2", d_b, 0)],
    );
    let (a, b) = (
        env.species.index_of("H_r").unwrap(),
        env.species.index_of("H2").unwrap(),
    );
    let (fine, coarse) = ([8, 2, 2], [4, 1, 1]);
    let (n0, runs, t_end) = (400i64, 400u64, 1.0e-6f64);

    // Mean field.
    let parent = |v: usize| {
        let c = [v % 8, (v / 8) % 2, v / 16];
        c[0] / 2
    };
    let rhs = |y: &[f64], dy: &mut [f64]| {
        dy.iter_mut().for_each(|x| *x = 0.0);
        let (ya, yb) = y.split_at(4);
        let (da, db) = dy.split_at_mut(4);
        add_laplacian(coarse, d_a / (2.0 * h).powi(2), ya, da);
        for (x, y) in da.iter_mut().zip(ya) {
            *x -= k * y;
        }
        add_laplacian(fine, d_b / (h * h), yb, db);
        for (v, x) in db.iter_mut().enumerate() {
            *x += k / 8.0 * ya[parent(v)];
        }
    };
    let mut y = vec![0.0; 36];
    y[0] = n0 as f64;
    let dt: f64 = 1e-10;
    let (mut k1, mut k2, mut k3, mut k4) =
        (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0.0; 36]);
    let mut tmp = vec![0.0; 36];
    for _ in 0..(t_end / dt).round() as usize {
        rhs(&y, &mut k1);
        tmp.iter_mut()
            .enumerate()
            .for_each(|(i, x)| *x = y[i] + 0.5 * dt * k1[i]);
        rhs(&tmp, &mut k2);
        tmp.iter_mut()
            .enumerate()
            .for_each(|(i, x)| *x = y[i] + 0.5 * dt * k2[i]);
        rhs(&tmp, &mut k3);
        tmp.iter_mut()
            .enumerate()
            .for_each(|(i, x)| *x = y[i] + dt * k3[i]);
        rhs(&tmp, &mut k4);
        for i in 0..36 {
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }

    // Ensemble.
    let (mut sum, mut sum_sq) = (vec![0.0; 36], vec![0.0; 36]);
    for stream in 0..runs {
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(11, stream));
        r.inject(&[(0, a, n0)]);
        r.advance_to(t_end);
        let t = r.totals();
        assert_eq!(t[a] + t[b], n0);
        let x: Vec<i64> = r.counts(a).iter().chain(r.counts(b)).copied().collect();
        for (i, n) in x.iter().enumerate() {
            sum[i] += *n as f64;
            sum_sq[i] += (*n as f64).powi(2);
        }
    }
    let n = runs as f64;
    let mut worst: f64 = 0.0;
    for i in 0..36 {
        let mean = sum[i] / n;
        let var = (sum_sq[i] - sum[i] * sum[i] / n) / (n - 1.0);
        let sem = (var / n).sqrt().max((y[i] / n).sqrt());
        let dev = (mean - y[i]).abs() / sem;
        worst = worst.max(dev);
        assert!(
            dev <= 4.0,
            "entry {i}: mean {mean} vs {} ({dev:.1} SEM)",
            y[i]
        );
    }
    // B has reached every fine voxel by then, and A has mostly converted.
    assert!(y[4..].iter().all(|x| *x > 0.05), "{:?}", &y[4..]);
    assert!(y[..4].iter().sum::<f64>() < 0.2 * n0 as f64);
    eprintln!("conversion: worst deviation {worst:.2} SEM");
}

#[test]
fn ssa_multigrid_cross_level_association_has_the_coarse_voxel_rate() {
    // One e_aq on the fine grid, one OH_r on level 2, which is a single
    // voxel: they react as if both occupied it (Hellander and Hellander
    // 2020, III.C.4, eq. 10-12), so the association time is exponential
    // with the (mesh-corrected) propensity of the coarse voxel, whatever
    // the hops of e_aq.
    let toml = "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = 5.0e10";
    let (h, d, k) = (2e-9, 5e-9, 5.0e10);
    let (env, s) = setup(
        toml,
        [4, 4, 4],
        h,
        2,
        true,
        &[("e_aq", d, 0), ("OH_r", d, 2), ("OH_minus", d, 0)],
    );
    let ix = |n: &str| env.species.index_of(n).unwrap();
    let (e, oh, ohm) = (ix("e_aq"), ix("OH_r"), ix("OH_minus"));

    // k [L/mol/s] -> m^3/s per pair; the level-2 voxel is 8 nm wide.
    let kappa = k / (1e3 * AVOGADRO);
    let h2 = 4.0 * h;
    let c_macro = kappa / h2.powi(3);
    let h_star = 1.5164 / 6.0 * kappa / (2.0 * d);
    let c_corr = c_macro / (1.0 - h_star / h2);
    assert!(c_corr > 1.3 * c_macro, "the correction must be visible");

    let runs = 4000u64;
    let (mut sum, mut sum_sq, mut hops) = (0.0, 0.0, 0u64);
    for stream in 0..runs {
        let mut r = Multigrid::new(&s, 0.0, SsaRng::new(21, stream));
        r.inject(&[((stream % 64) as usize, e, 1), (63, oh, 1)]);
        let mut at = (stream % 64) as usize;
        while r.n_reactions == 0 {
            at = r.counts(e).iter().position(|n| *n == 1).unwrap();
            assert!(r.step(f64::INFINITY));
        }
        // The product appears where the electron was.
        assert_eq!(r.counts(ohm)[at], 1);
        assert_eq!(r.totals()[e] + r.totals()[oh], 0);
        sum += r.t;
        sum_sq += r.t * r.t;
        hops += r.n_hops;
    }
    let n = runs as f64;
    let mean = sum / n;
    let sem = ((sum_sq - sum * sum / n) / (n - 1.0) / n).sqrt();
    let (t_corr, t_macro) = (1.0 / c_corr, 1.0 / c_macro);
    assert!(hops > 10 * runs, "e_aq hops: {hops}");
    assert!(
        (mean - t_corr).abs() <= 4.0 * sem,
        "mean {mean:.4e} s vs 1/k_meso {t_corr:.4e} s (SEM {sem:.1e})"
    );
    assert!((mean - t_macro).abs() > 10.0 * sem, "{mean} vs {t_macro}");
    eprintln!(
        "association: mean {mean:.4e} s, 1/k_meso {t_corr:.4e} s, {:.2} SEM",
        (mean - t_corr).abs() / sem
    );
}
