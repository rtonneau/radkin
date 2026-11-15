//! Reflective walls of the IRT stage: positions folded back into the box,
//! and the mirror images of a pair near a wall.

use std::path::{Path, PathBuf};

use radkin::irt::solver::fold;
use radkin::irt::{Boundary, IrtSimulation, Stage, StageOptions, run_irt_with};
use radkin::ssa::SsaRng;

/// Allowed distance in standard errors.
const K: f64 = 4.0;
const NM: f64 = 1e-9;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("radkin_irt_walls_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn toml_path(p: &Path) -> String {
    p.display().to_string().replace('\\', "/")
}

/// A run of `chemistry` in the reflective box `[0, edge_nm]^3`, with the
/// extra `[irt]` keys `irt` and the molecule CSV text `molecules`.
fn simulation(
    dir: &Path,
    name: &str,
    chemistry: &str,
    edge_nm: f64,
    irt: &str,
    molecules: &str,
) -> IrtSimulation {
    let mol = dir.join(format!("{name}.csv"));
    std::fs::write(&mol, molecules).unwrap();
    let text = format!(
        "{chemistry}\n[irt]\nboundary = \"reflective\"\n{irt}\n[irt.box]\nlo_nm = 0\nhi_nm = {edge_nm}\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
        toml_path(&mol)
    );
    let path = dir.join(format!("{name}.toml"));
    std::fs::write(&path, text).unwrap();
    IrtSimulation::from_toml(&path).unwrap()
}

#[test]
fn irt_walls_fold_reflects_into_the_box() {
    let (lo, hi) = ([0.0, -1.0, 2.0], [1.0, 1.0, 5.0]);
    assert_eq!(fold(lo, hi, [0.25, 0.5, 3.0]), [0.25, 0.5, 3.0]);
    // one reflection at the upper and at the lower wall
    assert_eq!(fold(lo, hi, [1.25, -1.5, 5.5]), [0.75, -0.5, 4.5]);
    assert_eq!(fold(lo, hi, [-0.25, 1.5, 1.0]), [0.25, 0.5, 3.0]);
    // several reflections
    assert_eq!(fold(lo, hi, [2.25, 3.5, 9.0]), [0.25, -0.5, 3.0]);
    assert_eq!(fold(lo, hi, [-3.75, -5.5, -5.0]), [0.25, -0.5, 3.0]);
    // on a wall
    assert_eq!(fold(lo, hi, [1.0, -1.0, 5.0]), [1.0, -1.0, 5.0]);
}

/// Chi-square of the positions over `slabs` slabs along each axis, against a
/// uniform distribution; returns the largest of the three.
fn slab_chi2(pos: &[[f64; 3]], edge: f64, slabs: usize) -> f64 {
    let expected = pos.len() as f64 / slabs as f64;
    (0..3)
        .map(|a| {
            let mut hits = vec![0usize; slabs];
            for p in pos {
                assert!((0.0..=edge).contains(&p[a]), "{p:?}");
                hits[((p[a] / edge * slabs as f64) as usize).min(slabs - 1)] += 1;
            }
            hits.iter()
                .map(|h| (*h as f64 - expected).powi(2) / expected)
                .sum::<f64>()
        })
        .fold(0.0, f64::max)
}

#[test]
fn irt_walls_uniform_distribution_stays_uniform() {
    let dir = scratch("uniform");
    // inert molecules: they only diffuse
    let chemistry = "[[reactions.k_reaction]]\nreaction = \"H2 -> H2O2\"\nk_value = 0";
    let sim = simulation(
        &dir,
        "inert",
        chemistry,
        100.0,
        "t_end = 1.0",
        "species,x,y,z\n",
    );
    let h2 = sim.model.species.iter().position(|s| s == "H2").unwrap() as u32;
    let edge = 100.0 * NM;
    let options = StageOptions {
        t_stop: 1.0,
        start_contact: false,
        boundary: Boundary::Reflective {
            lo: [0.0; 3],
            hi: [edge; 3],
        },
        cutoff_probability: 0.01,
        molecule_reactions: false,
        wall_images: true,
    };
    let n = 20_000;
    let slabs = 20;
    // chi-square with 19 degrees of freedom: mean 19, sd 6.2
    let limit = 19.0 + 5.0 * 38f64.sqrt();
    let mut rng = SsaRng::new(21, 0);
    let uniform: Vec<(u32, [f64; 3])> = (0..n)
        .map(|_| (h2, std::array::from_fn(|_| edge * rng.uniform())))
        .collect();
    let mut stage = Stage::new(&sim.model, options, SsaRng::new(22, 0), 0.0);
    stage.inject(&uniform);
    // moves from a tenth of the box to many boxes (D = 4.8e-9 m2/s)
    for t in [1e-9, 1e-8, 1e-7, 1e-6, 1e-4] {
        stage.sync(t);
        assert_eq!(stage.counts()[h2 as usize], n as i64);
        let chi2 = slab_chi2(&stage.particles().pos, edge, slabs);
        assert!(chi2 < limit, "t = {t:e} s: chi2 {chi2}");
    }
    // molecules that start in one corner fill the box
    let corner: Vec<(u32, [f64; 3])> = (0..n).map(|_| (h2, [0.1 * edge; 3])).collect();
    let mut stage = Stage::new(&sim.model, options, SsaRng::new(23, 0), 0.0);
    stage.inject(&corner);
    assert!(slab_chi2(&stage.particles().pos, edge, slabs) > 1e4);
    stage.sync(1e-5);
    let chi2 = slab_chi2(&stage.particles().pos, edge, slabs);
    assert!(chi2 < limit, "from a corner: chi2 {chi2}");
}

/// Survival of a pair in the box `[0, edge]^3` by a Brownian walk with
/// reflection: the fraction of `walks` pairs, started at `a` and `b`, that
/// have not come within `r_eff` by each time of `times`.
///
/// A step moves the pair by a tenth of its gap (of the reaction radius when
/// it is closer than that); a pair that crossed the radius between the two
/// ends of a step is caught by the crossing probability of a Brownian bridge.
fn brownian_survival(
    (a, b): ([f64; 3], [f64; 3]),
    (d_a, d_b): (f64, f64),
    r_eff: f64,
    edge: f64,
    times: &[f64],
    walks: usize,
) -> Vec<f64> {
    let d = d_a + d_b;
    let dist = |p: [f64; 3], q: [f64; 3]| {
        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    };
    let end = times[times.len() - 1];
    let mut alive = vec![0usize; times.len()];
    let (lo, hi) = ([0.0; 3], [edge; 3]);
    for w in 0..walks {
        let mut rng = SsaRng::new(31, w as u64);
        let (mut p, mut q) = (a, b);
        let mut t = 0.0;
        let mut d0 = dist(p, q);
        let mut reacted_at = f64::INFINITY;
        while t < end {
            let scale = (d0 - r_eff).max(r_eff);
            let dt = (0.005 * scale * scale / d).min(end - t);
            let (sa, sb) = ((2.0 * d_a * dt).sqrt(), (2.0 * d_b * dt).sqrt());
            for x in &mut p {
                *x += sa * rng.normal();
            }
            for x in &mut q {
                *x += sb * rng.normal();
            }
            p = fold(lo, hi, p);
            q = fold(lo, hi, q);
            t += dt;
            let d1 = dist(p, q);
            if d1 < r_eff || rng.uniform() < (-(d0 - r_eff) * (d1 - r_eff) / (d * dt)).exp() {
                reacted_at = t;
                break;
            }
            d0 = d1;
        }
        for (k, tk) in times.iter().enumerate() {
            if reacted_at > *tk {
                alive[k] += 1;
            }
        }
    }
    alive.iter().map(|n| *n as f64 / walks as f64).collect()
}

/// e_aq + OH_r, fully diffusion-controlled: r_eff 0.549 nm.
const PAIR_CHEMISTRY: &str = "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\n\
                              k_value = 2.95e10\nirt_type = \"full\"";

/// Survival of the e_aq of one e_aq + OH_r pair at `times`: the IRT stage
/// with and without the mirror images, and the Brownian walk.
struct PairSurvival {
    images: Vec<f64>,
    no_images: Vec<f64>,
    walk: Vec<f64>,
    trajectories: usize,
    walks: usize,
    image_reactions: u64,
    box_moves: u64,
}

fn pair_survival(
    name: &str,
    edge_nm: f64,
    (a, b): ([f64; 3], [f64; 3]),
    times: &[f64],
    trajectories: usize,
    walks: usize,
) -> PairSurvival {
    let dir = scratch(name);
    let molecules = format!(
        "species,x,y,z\ne_aq,{},{},{}\nOH_r,{},{},{}\n",
        a[0], a[1], a[2], b[0], b[1], b[2]
    );
    let list: Vec<String> = times.iter().map(|t| format!("{t:e}")).collect();
    let irt = format!(
        "t_end = {:e}\nstart_contact = false\nrecord_times = [{}]",
        1.01 * times[times.len() - 1],
        list.join(", ")
    );
    let run = |file: &str, irt: &str| {
        let sim = simulation(&dir, file, PAIR_CHEMISTRY, edge_nm, irt, &molecules);
        let out = run_irt_with(&sim, trajectories, 41).unwrap();
        let e = out.species.iter().position(|s| s == "e_aq").unwrap();
        // records: the start, the times, the end
        let s: Vec<f64> = (1..=times.len()).map(|k| out.totals_mean[k][e]).collect();
        (s, out.stats.clone(), sim)
    };
    let (images, stats, sim) = run("images", &irt);
    let (no_images, plain, _) = run("no_images", &format!("{irt}\nwall_images = false"));
    assert_eq!(plain.image_reactions, 0);
    let e = sim.model.species.iter().position(|s| s == "e_aq").unwrap();
    let oh = sim.model.species.iter().position(|s| s == "OH_r").unwrap();
    let walk = brownian_survival(
        (a.map(|x| x * NM), b.map(|x| x * NM)),
        (sim.model.diffusion[e], sim.model.diffusion[oh]),
        sim.model.pairs[0].r_eff,
        edge_nm * NM,
        times,
        walks,
    );
    println!("{name}: time [s]   walk     images   no images");
    for k in 0..times.len() {
        println!(
            "          {:.0e}     {:.4}   {:.4}   {:.4}",
            times[k], walk[k], images[k], no_images[k]
        );
    }
    PairSurvival {
        images,
        no_images,
        walk,
        trajectories,
        walks,
        image_reactions: stats.image_reactions,
        box_moves: stats.box_moves,
    }
}

fn binomial_sd(s: f64, n: usize) -> f64 {
    (s * (1.0 - s) / n as f64).sqrt()
}

/// A pair 4 nm from one wall of a 60 nm box (7 reaction radii from the
/// wall, the other walls out of reach): the mirror image of the partner is
/// 9.4 nm away, against 5 nm for the partner itself.
#[test]
fn irt_walls_pair_near_a_wall_against_a_brownian_walk() {
    let times = [1e-9, 5e-9, 2e-8];
    let s = pair_survival(
        "near_wall",
        60.0,
        ([4.0, 30.0, 30.0], [4.0, 30.0, 35.0]),
        &times,
        40_000,
        6_000,
    );
    assert!(s.image_reactions > 0);
    // one pass covers the run: the look-ahead limit of the box is not reached
    assert_eq!(s.box_moves, 0);
    for k in 0..times.len() {
        let err = (binomial_sd(s.walk[k], s.walks).powi(2)
            + binomial_sd(s.images[k], s.trajectories).powi(2))
        .sqrt();
        assert!(
            (s.images[k] - s.walk[k]).abs() <= K * err,
            "t = {:e} s: images {} against the walk {} (sd {err})",
            times[k],
            s.images[k],
            s.walk[k]
        );
    }
    // without the image the pair survives too often once the wall is felt
    let k = times.len() - 1;
    let err = (binomial_sd(s.walk[k], s.walks).powi(2)
        + binomial_sd(s.no_images[k], s.trajectories).powi(2))
    .sqrt();
    assert!(
        s.no_images[k] - s.walk[k] > K * err,
        "no images {} against the walk {} (sd {err})",
        s.no_images[k],
        s.walk[k]
    );
}

/// A pair in a box only 11 reaction radii wide: the limit of the image rule.
/// A sink and its image compete for the same molecules, which independent
/// pairs cannot see, so the rule reacts the pair too fast; without images it
/// reacts too slowly. The walk lies between the two.
#[test]
fn irt_walls_pair_in_a_tiny_box_is_bracketed() {
    let times = [3e-9, 6e-9, 1e-8];
    let s = pair_survival(
        "tiny_box",
        6.0,
        ([1.0, 1.5, 1.2], [4.6, 4.9, 4.4]),
        &times,
        4_000,
        2_000,
    );
    assert!(s.image_reactions > 0 && s.box_moves > 0);
    for k in 0..times.len() {
        assert!(
            s.images[k] < s.walk[k] && s.walk[k] < s.no_images[k],
            "t = {:e} s: images {}, walk {}, no images {}",
            times[k],
            s.images[k],
            s.walk[k],
            s.no_images[k]
        );
    }
    // the reaction rate of the image rule is within a factor 2 of the walk's
    let rate = |s: &[f64]| (s[0] / s[2]).ln() / (times[2] - times[0]);
    let ratio = rate(&s.images) / rate(&s.walk);
    assert!((1.0..2.0).contains(&ratio), "{ratio}");
}

/// Mean number of A left at `times` in a well-mixed volume, from `runs`
/// exact stochastic simulations of `A + B -> C` with `n0` of each and the
/// rate `c` per pair [1/s].
fn well_mixed(n0: u64, c: f64, times: &[f64], runs: usize) -> Vec<f64> {
    let mut sum = vec![0.0; times.len()];
    for r in 0..runs {
        let mut rng = SsaRng::new(51, r as u64);
        let mut n = n0;
        let mut t = 0.0;
        let mut k = 0;
        while k < times.len() {
            let next = if n == 0 {
                f64::INFINITY
            } else {
                t + rng.exp1() / (c * (n * n) as f64)
            };
            while k < times.len() && times[k] < next {
                sum[k] += n as f64;
                k += 1;
            }
            t = next;
            n = n.saturating_sub(1);
        }
    }
    sum.iter().map(|s| s / runs as f64).collect()
}

#[test]
fn irt_walls_well_mixed_box_matches_the_rate_equation() {
    let dir = scratch("mixed");
    // A + B -> C well below its diffusion limit (k_dif 1.5e10 M-1 s-1 for a
    // reaction radius of 0.5 nm and D = 2e-9 + 2e-9 m2/s)
    let k_obs = 3e9;
    let chemistry = "[[reactions.k_reaction]]\nreaction = \"O2 + O3 -> O_atom\"\nk_value = 3e9\n\
                     irt_type = \"partial\"\n\
                     [irt.species.O2]\nradius_nm = 0.25\ndiffusion = 2e-9\n\
                     [irt.species.O3]\nradius_nm = 0.25\ndiffusion = 2e-9";
    let edge_nm = 40.0;
    let n0 = 30;
    let mut rng = SsaRng::new(52, 0);
    let mut molecules = String::from("species,x,y,z\n");
    for i in 0..2 * n0 {
        let name = if i < n0 { "O2" } else { "O3" };
        let p: [f64; 3] = std::array::from_fn(|_| edge_nm * rng.uniform());
        molecules.push_str(&format!("{name},{},{},{}\n", p[0], p[1], p[2]));
    }
    let times = [2e-7, 5e-7, 1e-6, 1.5e-6];
    let irt = "t_end = 1.6e-6\nstart_contact = false\nrecord_times = [2e-7, 5e-7, 1e-6, 1.5e-6]";
    // the species table of the config goes with the chemistry text
    let (chem, table) = chemistry.split_once("[irt.species.O2]").unwrap();
    let irt = format!("{irt}\n[irt.species.O2]{table}");
    let sim = simulation(&dir, "mixed", chem, edge_nm, &irt, &molecules);
    let p = &sim.model.pairs[0];
    assert!((p.k_dif - 1.5e10).abs() < 0.02e10, "{}", p.k_dif);
    let trajectories = 100;
    let out = run_irt_with(&sim, trajectories, 53).unwrap();
    let a = out.species.iter().position(|s| s == "O2").unwrap();
    // rate per pair: k / (N_A V), V in litres
    let volume = (edge_nm * NM).powi(3) * 1e3;
    let c = k_obs / (6.02214076e23 * volume);
    let reference = well_mixed(n0 as u64, c, &times, 20_000);
    println!("time [s]   well mixed   irt");
    for k in 0..times.len() {
        let got = out.totals_mean[k + 1][a];
        let sem = out.totals_std[k + 1][a] / (trajectories as f64).sqrt();
        println!(
            "{:.1e}    {:.2}        {got:.2} +- {sem:.2}",
            times[k], reference[k]
        );
        assert!(
            (got - reference[k]).abs() <= K * sem + 0.02 * reference[k],
            "t = {:e} s: {got} against {} (sem {sem})",
            times[k],
            reference[k]
        );
    }
    // half of the molecules have reacted by the end
    assert!(reference[3] < 0.6 * n0 as f64 && reference[3] > 0.2 * n0 as f64);
}

#[test]
fn irt_walls_seed_outside_the_box_is_refused() {
    let dir = scratch("seed");
    let chemistry =
        "[[reactions.k_reaction]]\nreaction = \"e_aq + OH_r -> OH_minus\"\nk_value = 2.95e10";
    let mol = dir.join("m.csv");
    std::fs::write(&mol, "species,x,y,z\ne_aq,1,1,1\nOH_r,7,1,1\n").unwrap();
    let text = format!(
        "{chemistry}\n[irt]\nt_end = 1e-9\nboundary = \"reflective\"\n[irt.box]\nlo_nm = 0\nhi_nm = 6\n\
         [[irt.seed]]\nfile = \"{}\"\ntime = 0\n",
        toml_path(&mol)
    );
    let path = dir.join("bad.toml");
    std::fs::write(&path, &text).unwrap();
    let e = IrtSimulation::from_toml(&path).unwrap_err().to_string();
    assert!(e.contains("outside the box"), "{e}");
    // without a box
    let no_box = text.replace("[irt.box]\nlo_nm = 0\nhi_nm = 6\n", "");
    std::fs::write(&path, no_box).unwrap();
    let e = IrtSimulation::from_toml(&path).unwrap_err().to_string();
    assert!(e.contains("needs a box"), "{e}");
}
