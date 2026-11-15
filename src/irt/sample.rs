//! Sampling kernel of the IRT stage: the reaction time of one pair, the
//! reaction of a pair found in contact, and the positions of the reactants
//! and products of a reaction.
//!
//! The formulas are those of Geant4-DNA 11.4.3
//! (`G4DiffusionControlledReactionModel::GetTimeToEncounter`,
//! `G4DNAMakeReaction`). Two points differ on purpose, where Geant4's code
//! departs from the formulas it implements (Plante and Devroye 2017):
//!
//! - the diffusion of the pair during the step uses the standard deviation
//!   `sqrt(sigma)`; Geant4 passes the variance `sigma` to its Gaussian;
//! - the angle of the separation at the reaction is taken from the initial
//!   separation; Geant4 takes it from the z axis.
//!
//! The kernel does not allocate and uses IEEE basic operations, the square
//! root, [`det_ln`] and [`det_exp`] only, so a seed gives the same draws on
//! every platform.

use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

use crate::chemistry::reactions::IrtType;
use crate::irt::model::PairReaction;
use crate::ssa::rng::{SsaRng, det_exp, det_ln, norm_quantile};

/// A position or a displacement [m].
pub type Vec3 = [f64; 3];

/// Inverse of the complementary error function, for `x` in `(0, 2)`.
pub fn erfc_inv(x: f64) -> f64 {
    -FRAC_1_SQRT_2 * norm_quantile(0.5 * x)
}

/// Reaction time [s] of a pair `distance` apart, from now; `None` when the
/// pair never reacts.
///
/// The pair reaches the effective radius with the probability
/// `W = r_eff / distance`, at the time `(distance - r_eff)^2 / (4 D
/// erfcinv(U / W)^2)`. A partially diffusion-controlled pair then reacts with
/// the probability `k_act / (k_act + k_dif)`, after a second time. A
/// `distance` below `r_eff` counts as a pair in contact.
pub fn pair_time(rng: &mut SsaRng, p: &PairReaction, distance: f64) -> Option<f64> {
    let w = p.r_eff / distance;
    let u = rng.uniform_inside();
    if u >= w {
        return None;
    }
    let t1 = if w >= 1.0 {
        0.0
    } else {
        let e = erfc_inv(u / w);
        if e <= 0.0 {
            return None;
        }
        let d = (distance - p.r_eff) / e;
        d * d / (4.0 * p.d_sum)
    };
    if p.kind == IrtType::Full {
        return Some(t1);
    }
    let sum = p.k_act + p.k_dif;
    if rng.uniform() > p.k_act / sum {
        return None;
    }
    let y = (SQRT_2 * rng.normal()).abs();
    let x = if y > 0.0 { rng.exp1() / y } else { 0.0 };
    let f = x * p.r_eff * p.k_dif / sum;
    Some(t1 + f * f / p.d_sum)
}

/// Whether a pair found inside the effective radius at the start reacts at
/// once.
pub fn contact_reacts(rng: &mut SsaRng, p: &PairReaction) -> bool {
    p.p_contact >= 1.0 || rng.uniform() < p.p_contact
}

fn norm(v: Vec3) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Three independent standard normal draws.
fn normal3(rng: &mut SsaRng) -> Vec3 {
    [rng.normal(), rng.normal(), rng.normal()]
}

/// Uniform direction.
fn direction(rng: &mut SsaRng) -> Vec3 {
    loop {
        let g = normal3(rng);
        let n = norm(g);
        if n > 1e-12 {
            return [g[0] / n, g[1] / n, g[2] / n];
        }
    }
}

/// Cosine of the angle between the separation of a pair at its reaction and
/// its initial separation, for `alpha = r_eff * distance / (2 D dt)`.
fn cos_theta(rng: &mut SsaRng, alpha: f64) -> f64 {
    let u = rng.uniform();
    let c = if alpha < 1e-8 {
        1.0 - 2.0 * u
    } else {
        1.0 + det_ln(1.0 - u * (1.0 - det_exp(-2.0 * alpha))) / alpha
    };
    c.clamp(-1.0, 1.0)
}

/// Positions of the two reactants at the reaction, `dt` after they were at
/// `pa` and `pb`.
///
/// The separation has the length `r_eff`. Its direction makes an angle with
/// the initial separation drawn for `dt`, and the coordinate
/// `pa + (D_a / D_b) pb`, which diffuses independently of the separation,
/// moves by its free diffusion over `dt`. A reactant that does not diffuse
/// stays where it is and the other one joins it.
pub fn reaction_positions(
    rng: &mut SsaRng,
    r_eff: f64,
    (d_a, d_b): (f64, f64),
    (pa, pb): (Vec3, Vec3),
    dt: f64,
) -> (Vec3, Vec3) {
    if d_a == 0.0 {
        return (pa, pa);
    }
    if d_b == 0.0 {
        return (pb, pb);
    }
    if dt <= 0.0 {
        return (pa, pb);
    }
    let sep = [pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]];
    let distance = norm(sep);
    let d = d_a + d_b;
    let axis = if distance > 0.0 {
        [sep[0] / distance, sep[1] / distance, sep[2] / distance]
    } else {
        direction(rng)
    };
    let cos = cos_theta(rng, r_eff * distance / (2.0 * d * dt));
    let sin = (1.0 - cos * cos).max(0.0).sqrt();
    // Unit vector perpendicular to the axis, uniform around it.
    let perp = loop {
        let g = normal3(rng);
        let along = g[0] * axis[0] + g[1] * axis[1] + g[2] * axis[2];
        let v = [
            g[0] - along * axis[0],
            g[1] - along * axis[1],
            g[2] - along * axis[2],
        ];
        let n = norm(v);
        if n > 1e-12 {
            break [v[0] / n, v[1] / n, v[2] / n];
        }
    };
    let s12 = 2.0 * d_a * dt;
    let s22 = 2.0 * d_b * dt;
    let sd = (s12 + s12 * s12 / s22).sqrt();
    let ratio = d_a / d_b;
    let mut ra = [0.0; 3];
    let mut rb = [0.0; 3];
    for i in 0..3 {
        let s1 = r_eff * (cos * axis[i] + sin * perp[i]);
        let s2 = pa[i] + ratio * pb[i] + sd * rng.normal();
        ra[i] = (d_a * s1 + d_b * s2) / d;
        rb[i] = d_b * (s2 - s1) / d;
    }
    (ra, rb)
}

/// Position of the products of a reaction between reactants at `pa` and
/// `pb`: the mean of the diffusion-weighted reaction site and a uniform point
/// of the segment between the reactants. Every product of one reaction gets
/// the same position.
pub fn product_position(rng: &mut SsaRng, (d_a, d_b): (f64, f64), (pa, pb): (Vec3, Vec3)) -> Vec3 {
    let (sa, sb) = (d_a.sqrt(), d_b.sqrt());
    let inv = 1.0 / (sa + sb);
    let u = rng.uniform();
    let mut out = [0.0; 3];
    for i in 0..3 {
        let site = (sb * pa[i] + sa * pb[i]) * inv;
        let random = (1.0 - u) * pa[i] + u * pb[i];
        out[i] = 0.5 * (site + random);
    }
    out
}

/// Free diffusion of a molecule over `dt`.
pub fn diffuse(rng: &mut SsaRng, d: f64, dt: f64, pos: &mut Vec3) {
    if d <= 0.0 || dt <= 0.0 {
        return;
    }
    let sd = (2.0 * d * dt).sqrt();
    for x in pos {
        *x += sd * rng.normal();
    }
}
