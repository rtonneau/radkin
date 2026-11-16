//! Sampling kernel of the IRT stage against analytic survival probabilities.

use radkin::chemistry::reactions::IrtType;
use radkin::irt::model::{PairReaction, TypeSource, effective_radius, onsager_radius};
use radkin::irt::sample::{
    contact_reacts, diffuse, erfc_inv, pair_time, product_position, reaction_positions,
};
use radkin::ssa::SsaRng;

const AVOGADRO: f64 = 6.02214076e23;
const PI: f64 = std::f64::consts::PI;

/// Allowed distance in standard errors.
const K: f64 = 4.0;

/// Complementary error function (reference of the tests): the series of
/// `erf` below 2, a continued fraction above.
fn erfc(z: f64) -> f64 {
    if z < 0.0 {
        return 2.0 - erfc(-z);
    }
    if z < 2.0 {
        // erf z = 2/sqrt(pi) exp(-z^2) sum_n 2^n z^(2n+1) / (1 3 ... (2n+1))
        let mut term = z;
        let mut sum = z;
        let mut n = 0.0;
        while term > 1e-18 * sum {
            n += 1.0;
            term *= 2.0 * z * z / (2.0 * n + 1.0);
            sum += term;
        }
        return 1.0 - 2.0 / PI.sqrt() * (-z * z).exp() * sum;
    }
    (-z * z).exp() * erfcx(z)
}

/// `exp(z^2) erfc(z)` for `z >= 0`, without overflow.
fn erfcx(z: f64) -> f64 {
    if z < 2.0 {
        return (z * z).exp() * erfc(z);
    }
    // 1 / (z sqrt(pi)) / (1 + 1/(2z^2) / (1 + 2/(2z^2) / (1 + ...)))
    let a = 1.0 / (2.0 * z * z);
    let mut f = 1.0;
    for n in (1..=120).rev() {
        f = 1.0 + n as f64 * a / f;
    }
    1.0 / (z * PI.sqrt()) / f
}

/// Pair reaction of two species with the diffusion coefficients `d` [m^2/s].
fn full(k_obs: f64, d: (f64, f64)) -> PairReaction {
    let d_sum = d.0 + d.1;
    let r_eff = k_obs / 1e3 / (4.0 * PI * d_sum * AVOGADRO);
    PairReaction {
        index: 0,
        a: 0,
        b: 1,
        kind: IrtType::Full,
        source: TypeSource::Network,
        k_obs,
        d_sum,
        r_reaction: r_eff,
        r_onsager: 0.0,
        r_eff,
        k_dif: k_obs,
        k_act: f64::INFINITY,
        p_contact: 1.0,
        products: Vec::new(),
    }
}

fn partial(k_obs: f64, d: (f64, f64), radius_sum: f64, charges: (i32, i32)) -> PairReaction {
    let d_sum = d.0 + d.1;
    let rc = onsager_radius(charges.0, charges.1);
    let r_eff = effective_radius(radius_sum, rc);
    let k_dif = 4.0 * PI * d_sum * r_eff * AVOGADRO * 1e3;
    assert!(k_obs < k_dif);
    let k_act = k_dif * k_obs / (k_dif - k_obs);
    PairReaction {
        index: 0,
        a: 0,
        b: 1,
        kind: IrtType::Partial,
        source: TypeSource::Network,
        k_obs,
        d_sum,
        r_reaction: radius_sum,
        r_onsager: rc,
        r_eff,
        k_dif,
        k_act,
        p_contact: 0.29e-9 / (0.29e-9 + (k_dif / k_act) * (r_eff + 0.29e-9)),
        products: Vec::new(),
    }
}

/// Probability that a pair `r` apart has reacted by `t`.
fn reacted(p: &PairReaction, r: f64, t: f64) -> f64 {
    let rho = r - p.r_eff;
    let x = rho / (4.0 * p.d_sum * t).sqrt();
    if p.kind == IrtType::Full {
        return p.r_eff / r * erfc(x);
    }
    // radiation boundary: beta = (k_act + k_dif) / (k_dif r_eff)
    let beta = (p.k_act + p.k_dif) / (p.k_dif * p.r_eff);
    let y = beta * (p.d_sum * t).sqrt();
    p.r_eff / r
        * p.k_act
        / (p.k_act + p.k_dif)
        // exp(beta rho + y^2) erfc(x + y) = exp(-x^2) erfcx(x + y), as 2 x y = beta rho
        * (erfc(x) - (-x * x).exp() * erfcx(x + y))
}

/// Sampled reaction times of `n` pairs `r` apart (`None`: no reaction).
fn sample(p: &PairReaction, r: f64, n: usize, seed: u64) -> Vec<Option<f64>> {
    let mut rng = SsaRng::new(seed, 0);
    (0..n).map(|_| pair_time(&mut rng, p, r)).collect()
}

fn assert_fraction(hits: usize, n: usize, expected: f64, what: &str) {
    let got = hits as f64 / n as f64;
    let se = (expected * (1.0 - expected) / n as f64).sqrt();
    assert!(
        (got - expected).abs() <= K * se,
        "{what}: {got} against {expected} (se {se})"
    );
}

/// Empirical reaction probability at several times against `reacted`.
fn assert_survival(p: &PairReaction, r: f64, times: &[f64], seed: u64, what: &str) {
    let n = 400_000;
    let draws = sample(p, r, n, seed);
    for &t in times {
        let hits = draws.iter().filter(|x| x.is_some_and(|x| x <= t)).count();
        assert_fraction(hits, n, reacted(p, r, t), &format!("{what}, t = {t:e} s"));
    }
    // long-time limit
    let total = draws.iter().filter(|x| x.is_some()).count();
    let limit = match p.kind {
        IrtType::Full => p.r_eff / r,
        IrtType::Partial => p.r_eff / r * p.k_act / (p.k_act + p.k_dif),
    };
    assert_fraction(total, n, limit, &format!("{what}, long time"));
}

#[test]
fn irt_sample_erfc_inv_inverts_erfc() {
    let mut x = 1e-300;
    while x < 1.0 {
        let z = erfc_inv(x);
        let back = erfc(z);
        assert!(
            (back - x).abs() <= 1e-12 * x,
            "erfc(erfc_inv({x})) = {back}"
        );
        x *= 3.7;
    }
    for i in 1..400 {
        let x = i as f64 / 200.0;
        let back = erfc(erfc_inv(x));
        assert!((back - x).abs() <= 1e-12, "erfc(erfc_inv({x})) = {back}");
    }
    assert_eq!(erfc_inv(1.0), 0.0);
    assert!((erfc_inv(0.5) - 0.476_936_276_204_469_9).abs() < 1e-15);
    assert!((erfc_inv(1.5) + 0.476_936_276_204_469_9).abs() < 1e-15);
}

#[test]
fn irt_sample_fully_diffusion_controlled_survival() {
    // e_aq + OH_r like: k 2.95e10, D 4.9e-9 + 2.2e-9
    let p = full(2.95e10, (4.9e-9, 2.2e-9));
    assert!((p.r_eff - 0.549e-9).abs() < 1e-12, "{}", p.r_eff);
    for (r, seed) in [(0.8e-9, 1), (2.0e-9, 2), (10.0e-9, 3)] {
        assert_survival(
            &p,
            r,
            &[1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6],
            seed,
            &format!("full, r = {r:e}"),
        );
    }
}

#[test]
fn irt_sample_partially_diffusion_controlled_survival() {
    // OH_r + OH_r like: k_obs 5.5e9 on distinct species, radii 0.22 + 0.22 nm
    let p = partial(5.5e9, (2.2e-9, 2.2e-9), 0.44e-9, (0, 0));
    assert_eq!(p.r_eff, 0.44e-9);
    for (r, seed) in [(0.6e-9, 4), (1.5e-9, 5), (6.0e-9, 6)] {
        assert_survival(
            &p,
            r,
            &[1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6],
            seed,
            &format!("partial, r = {r:e}"),
        );
    }
}

#[test]
fn irt_sample_ion_pair_survival() {
    // e_aq + H3O_plus: k_obs 2.11e10, radii 0.5 + 0.25 nm, opposite charges
    let p = partial(2.11e10, (4.9e-9, 9.46e-9), 0.75e-9, (-1, 1));
    assert!((p.r_eff - 1.161_259_398e-9).abs() < 1e-18, "{}", p.r_eff);
    assert!(p.r_onsager < 0.0);
    for (r, seed) in [(1.5e-9, 7), (5.0e-9, 8)] {
        assert_survival(
            &p,
            r,
            &[1e-12, 1e-11, 1e-10, 1e-9, 1e-7],
            seed,
            &format!("attractive ions, r = {r:e}"),
        );
    }
    // like charges: the effective radius is below the reaction radius
    let p = partial(1e8, (2.0e-9, 2.0e-9), 0.5e-9, (-1, -1));
    assert!(p.r_eff < 0.5e-9 && p.r_onsager > 0.0);
    assert_survival(&p, 1.0e-9, &[1e-10, 1e-8], 9, "repulsive ions");
}

#[test]
fn irt_sample_pair_in_contact() {
    // a pair at the effective radius reaches it at once
    let p = full(2.95e10, (4.9e-9, 2.2e-9));
    let mut rng = SsaRng::new(10, 0);
    for _ in 0..100 {
        assert_eq!(pair_time(&mut rng, &p, p.r_eff), Some(0.0));
        assert_eq!(pair_time(&mut rng, &p, 0.5 * p.r_eff), Some(0.0));
        assert!(contact_reacts(&mut rng, &p));
    }
    // contact probability of a partially diffusion-controlled pair
    let p = partial(5.5e9, (2.2e-9, 2.2e-9), 0.44e-9, (0, 0));
    let n = 200_000;
    let hits = (0..n).filter(|_| contact_reacts(&mut rng, &p)).count();
    assert_fraction(hits, n, p.p_contact, "contact probability");
    // in contact, it reacts later with the probability k_act / (k_act + k_dif)
    let draws: Vec<_> = (0..n).map(|_| pair_time(&mut rng, &p, p.r_eff)).collect();
    let hits = draws.iter().filter(|x| x.is_some()).count();
    assert_fraction(hits, n, p.k_act / (p.k_act + p.k_dif), "contact, long time");
    assert!(draws.iter().flatten().all(|&t| t >= 0.0));
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[test]
fn irt_sample_reactants_meet_at_the_effective_radius() {
    let mut rng = SsaRng::new(12, 0);
    let (da, db): (f64, f64) = (4.9e-9, 2.2e-9);
    let r_eff = 0.55e-9;
    let pa = [1.0e-9, 2.0e-9, -3.0e-9];
    let pb = [2.5e-9, 2.0e-9, -1.0e-9];
    let dt: f64 = 2e-10;
    let n = 100_000;
    let ratio = da / db;
    // variance per axis of pa + (D_a / D_b) pb over dt
    let var = 2.0 * da * dt + (2.0 * da * dt).powi(2) / (2.0 * db * dt);
    let (mut m1, mut m2, mut cos_sum) = ([0.0; 3], [0.0; 3], 0.0);
    let d0 = dist(pa, pb);
    for _ in 0..n {
        let (ra, rb) = reaction_positions(&mut rng, r_eff, (da, db), (pa, pb), dt);
        assert!((dist(ra, rb) - r_eff).abs() < 1e-9 * r_eff);
        for i in 0..3 {
            let w = ra[i] + ratio * rb[i] - (pa[i] + ratio * pb[i]);
            m1[i] += w;
            m2[i] += w * w;
            cos_sum += (ra[i] - rb[i]) * (pa[i] - pb[i]) / (r_eff * d0);
        }
    }
    let nf = n as f64;
    for i in 0..3 {
        let mean = m1[i] / nf;
        let v = m2[i] / nf;
        assert!(mean.abs() < K * (var / nf).sqrt(), "axis {i}: mean {mean}");
        assert!(
            (v - var).abs() < K * var * (2.0 / nf).sqrt(),
            "axis {i}: {v}"
        );
    }
    // mean cosine of the angle to the initial separation: coth(a) - 1/a
    let alpha = r_eff * d0 / (2.0 * (da + db) * dt);
    let expected = 1.0 / alpha.tanh() - 1.0 / alpha;
    let got = cos_sum / nf;
    assert!(
        (got - expected).abs() < K / nf.sqrt(),
        "cos {got} against {expected}"
    );
}

#[test]
fn irt_sample_reaction_positions_special_cases() {
    let mut rng = SsaRng::new(13, 0);
    let pa = [1.0e-9, 0.0, 0.0];
    let pb = [0.0, 2.0e-9, 0.0];
    // contact reaction: nobody moves
    assert_eq!(
        reaction_positions(&mut rng, 0.5e-9, (1e-9, 2e-9), (pa, pb), 0.0),
        (pa, pb)
    );
    // a reactant that does not diffuse is joined by the other one
    assert_eq!(
        reaction_positions(&mut rng, 0.5e-9, (0.0, 2e-9), (pa, pb), 1e-9),
        (pa, pa)
    );
    assert_eq!(
        reaction_positions(&mut rng, 0.5e-9, (1e-9, 0.0), (pa, pb), 1e-9),
        (pb, pb)
    );
    // a long step forgets the initial direction: the mean cosine goes to 0
    let n = 50_000;
    let mut cos_sum = 0.0;
    for _ in 0..n {
        let (ra, rb) = reaction_positions(&mut rng, 0.5e-9, (1e-9, 2e-9), (pa, pb), 1.0);
        for i in 0..3 {
            cos_sum += (ra[i] - rb[i]) * (pa[i] - pb[i]) / (0.5e-9 * dist(pa, pb));
        }
    }
    assert!((cos_sum / n as f64).abs() < K / (3.0 * n as f64).sqrt());
}

#[test]
fn irt_sample_products_between_the_reactants() {
    let mut rng = SsaRng::new(14, 0);
    let pa = [0.0, 0.0, 0.0];
    let pb = [1.0e-9, 0.0, 0.0];
    let (da, db) = (4.0e-9, 1.0e-9);
    // site: (sqrt(db) pa + sqrt(da) pb) / (sqrt(da) + sqrt(db)) = 2/3 nm
    let n = 100_000;
    let (mut lo, mut hi, mut sum) = (f64::MAX, f64::MIN, 0.0);
    for _ in 0..n {
        let p = product_position(&mut rng, (da, db), (pa, pb));
        assert_eq!((p[1], p[2]), (0.0, 0.0));
        lo = lo.min(p[0]);
        hi = hi.max(p[0]);
        sum += p[0];
    }
    // x = (2/3 + u) / 2 nm with u uniform in [0, 1)
    assert!((1.0e-9 / 3.0..0.334e-9).contains(&lo), "{lo}");
    assert!((0.833e-9..2.5e-9 / 3.0).contains(&hi), "{hi}");
    let mean = sum / n as f64;
    assert!((mean - 7.0e-9 / 12.0).abs() < K * 0.5e-9 / (12.0 * n as f64).sqrt());
}

#[test]
fn irt_sample_free_diffusion_variance() {
    let mut rng = SsaRng::new(15, 0);
    let (d, dt) = (2.2e-9, 1e-9);
    let n = 100_000;
    let mut m2 = [0.0; 3];
    for _ in 0..n {
        let mut p = [0.0; 3];
        diffuse(&mut rng, d, dt, &mut p);
        for i in 0..3 {
            m2[i] += p[i] * p[i];
        }
    }
    let var = 2.0 * d * dt;
    for v in m2 {
        let v = v / n as f64;
        assert!((v - var).abs() < K * var * (2.0 / n as f64).sqrt());
    }
    let mut p = [1.0, 2.0, 3.0];
    diffuse(&mut rng, 0.0, dt, &mut p);
    diffuse(&mut rng, d, 0.0, &mut p);
    assert_eq!(p, [1.0, 2.0, 3.0]);
}
