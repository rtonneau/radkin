//! Event loop of the IRT stage with incremental resampling.
//!
//! Every sampled reaction time is kept in a queue. A reaction time is drawn
//! for a pair once, when the later of its two particles appears (injected,
//! or created by a reaction), and for a particle alone when it appears. When
//! a reaction fires, only the pairs of its products are drawn; the times of
//! every other pair stay as they are.
//!
//! A particle keeps the position of the last time it was looked at. It is
//! moved to the present, by its free diffusion since then, when a new
//! particle looks for partners around it, and when it reacts alone. A pair
//! that fires is put at its reaction positions from the two positions its
//! time was drawn from.
//!
//! The rules shared with the full resampling ([`crate::irt::solver`]): the
//! contact rule for just injected particles, a pair found inside its
//! effective radius later is looked at again after `r_eff^2 / D`, the
//! outflow box, and the mirror images between reflective walls. Between
//! reflective walls every pair is drawn again at the look-ahead limit of the
//! box, as the full resampling does at its box moves: only then can two
//! particles that started far apart find each other.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::irt::model::IrtModel;
use crate::irt::sample::{
    Vec3, contact_reacts, diffuse, erfc_inv, pair_time, product_position, reaction_positions,
};
use crate::irt::search::CellGrid;
use crate::irt::solver::{Boundary, StageOptions, StageStats, fold, mirror, outside};
use crate::ssa::rng::SsaRng;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Kind {
    /// Reaction of the particles `i` and `j`, drawn at `t_sample` from the
    /// positions `at`.
    Pair {
        i: u32,
        j: u32,
        reaction: u32,
        image: [i8; 3],
        t_sample: f64,
        at: (Vec3, Vec3),
        /// A pair in contact at its injection.
        contact: bool,
        epoch: u32,
    },
    /// Reaction of the particle `i` alone.
    Molecule { i: u32, reaction: u32 },
    /// Look again at a pair found inside its effective radius.
    Relook { i: u32, j: u32, epoch: u32 },
    /// Draw every pair again (reflective walls).
    BoxPass,
}

#[derive(Debug, Clone, Copy)]
struct Event {
    time: f64,
    /// Order of creation: breaks the ties.
    seq: u64,
    kind: Kind,
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Event {}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    /// The earliest event is the greatest: `BinaryHeap` is a max-heap.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .time
            .total_cmp(&self.time)
            .then(other.seq.cmp(&self.seq))
    }
}

/// One trajectory of the IRT stage with incremental resampling.
#[derive(Debug, Clone)]
pub struct IncrementalStage<'a> {
    model: &'a IrtModel,
    opt: StageOptions,
    rng: SsaRng,
    /// Every particle ever added; a particle keeps its index.
    species: Vec<u32>,
    pos: Vec<Vec3>,
    /// Time of each position [s].
    t_pos: Vec<f64>,
    alive: Vec<bool>,
    counts: Vec<i64>,
    queue: BinaryHeap<Event>,
    seq: u64,
    /// Pair events of an earlier epoch are void.
    epoch: u32,
    /// Stage time [s]: the last event fired, or the last `sync`.
    t: f64,
    box_pass_scheduled: bool,
    stats: StageStats,
    grid: CellGrid,
}

fn distance(a: Vec3, b: Vec3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

impl<'a> IncrementalStage<'a> {
    /// Empty stage at `t_start`.
    pub fn new(model: &'a IrtModel, opt: StageOptions, rng: SsaRng, t_start: f64) -> Self {
        let n = model.n_species();
        Self {
            model,
            opt,
            rng,
            species: Vec::new(),
            pos: Vec::new(),
            t_pos: Vec::new(),
            alive: Vec::new(),
            counts: vec![0; n],
            queue: BinaryHeap::new(),
            seq: 0,
            epoch: 0,
            t: t_start,
            box_pass_scheduled: false,
            stats: StageStats {
                by_reaction: vec![0; model.n_reactions()],
                absorbed: vec![0; n],
                background_consumed: vec![0; n],
                ..StageStats::default()
            },
            grid: CellGrid::new(),
        }
    }

    /// Stage time [s].
    pub fn time(&self) -> f64 {
        self.t
    }

    /// Particles per species.
    pub fn counts(&self) -> &[i64] {
        &self.counts
    }

    pub fn stats(&self) -> &StageStats {
        &self.stats
    }

    /// `(species, position)` of every particle, each at the time it was last
    /// looked at ([`IncrementalStage::sync`] brings them all to one time).
    pub fn molecules(&self) -> Vec<(u32, Vec3)> {
        (0..self.species.len())
            .filter(|i| self.alive[*i])
            .map(|i| (self.species[i], self.pos[i]))
            .collect()
    }

    fn push(&mut self, time: f64, kind: Kind) {
        self.seq += 1;
        self.queue.push(Event {
            time,
            seq: self.seq,
            kind,
        });
    }

    fn add(&mut self, species: u32, pos: Vec3, t: f64) -> usize {
        self.species.push(species);
        self.pos.push(pos);
        self.t_pos.push(t);
        self.alive.push(true);
        self.counts[species as usize] += 1;
        self.species.len() - 1
    }

    fn kill(&mut self, i: usize) {
        self.alive[i] = false;
        self.counts[self.species[i] as usize] -= 1;
    }

    /// Move the particle `i` to `t` by its free diffusion; false if it left
    /// the box on the way.
    fn move_to(&mut self, i: usize, t: f64) -> bool {
        let dt = t - self.t_pos[i];
        if dt > 0.0 {
            let s = self.species[i] as usize;
            diffuse(&mut self.rng, self.model.diffusion[s], dt, &mut self.pos[i]);
            self.t_pos[i] = t;
            match self.opt.boundary {
                Boundary::Outflow { lo, hi } if outside(lo, hi, self.pos[i]) => {
                    self.kill(i);
                    self.stats.absorbed[s] += 1;
                    return false;
                }
                Boundary::Reflective { lo, hi } => self.pos[i] = fold(lo, hi, self.pos[i]),
                _ => {}
            }
        }
        true
    }

    /// Look-ahead of a pair drawn at `t` [s], and the look-ahead limit of a
    /// box with reflective walls if there is one.
    fn horizon(&self, t: f64) -> (f64, Option<f64>) {
        let left = self.opt.t_stop - t;
        match self.opt.boundary {
            Boundary::Reflective { lo, hi } if self.model.max_d_sum() > 0.0 => {
                let edge = (0..3).map(|a| hi[a] - lo[a]).fold(f64::INFINITY, f64::min);
                let reach = edge / (2.0 * erfc_inv(self.opt.cutoff_probability));
                let limit = reach * reach / self.model.max_d_sum();
                (left.min(limit), Some(limit))
            }
            _ => (left, None),
        }
    }

    fn cutoff(&self, horizon: f64) -> f64 {
        self.model.max_radius()
            + 2.0
                * (self.model.max_d_sum() * horizon).sqrt()
                * erfc_inv(self.opt.cutoff_probability)
    }

    /// Add `(species, position)` particles at the stage time and draw their
    /// reaction times.
    pub fn inject(&mut self, molecules: &[(u32, Vec3)]) {
        let t = self.t;
        let first = self.species.len();
        for &(s, pos) in molecules {
            self.add(s, pos, t);
        }
        for p in first..self.species.len() {
            self.draw(p, t, self.opt.start_contact);
        }
        if let (false, (_, Some(limit))) = (self.box_pass_scheduled, self.horizon(t))
            && t + limit < self.opt.t_stop
        {
            self.box_pass_scheduled = true;
            self.push(t + limit, Kind::BoxPass);
        }
    }

    /// Fire every reaction up to `t`. The particles stay where they were
    /// last looked at: only their counts are those of `t`.
    pub fn advance_to(&mut self, t: f64) {
        while let Some(top) = self.queue.peek() {
            if top.time > t {
                break;
            }
            let event = self.queue.pop().expect("an event");
            self.fire(event);
        }
    }

    /// Move every particle to `t` by its free diffusion. Reactions up to `t`
    /// must have fired ([`IncrementalStage::advance_to`]). The reaction times
    /// already drawn are kept.
    pub fn sync(&mut self, t: f64) {
        for i in 0..self.species.len() {
            if self.alive[i] {
                self.move_to(i, t);
            }
        }
        self.t = self.t.max(t);
    }

    /// Draw the reaction times of the particle `p`, which appeared at `t`:
    /// alone, and with every particle that was there before it.
    fn draw(&mut self, p: usize, t: f64, fresh: bool) {
        let model = self.model;
        let sp = self.species[p] as usize;
        if self.opt.molecule_reactions {
            let rate = model.molecule_rate(sp);
            if rate > 0.0 {
                let dt = self.rng.exp1() / rate;
                if t + dt < self.opt.t_stop {
                    let list = model.molecule_of(sp);
                    let mut v = self.rng.uniform() * rate;
                    let mut reaction = list[list.len() - 1];
                    for &r in list {
                        let rate = model.molecule[r as usize].rate;
                        if v < rate {
                            reaction = r;
                            break;
                        }
                        v -= rate;
                    }
                    self.push(
                        t + dt,
                        Kind::Molecule {
                            i: p as u32,
                            reaction,
                        },
                    );
                }
            }
        }
        let (horizon, _) = self.horizon(t);
        if horizon.is_nan() || horizon <= 0.0 {
            return;
        }
        let cutoff = self.cutoff(horizon);
        for q in 0..p {
            if !self.alive[q] || !self.alive[p] {
                continue;
            }
            let sq = self.species[q] as usize;
            if model.pairs_of(sp, sq).is_empty() {
                continue;
            }
            // `q` may be where it was some time ago: it can have come closer
            // by 5 standard deviations of its displacement since then.
            let age = t - self.t_pos[q];
            let slack = if age > 0.0 {
                5.0 * (2.0 * model.diffusion[sq] * age).sqrt()
            } else {
                0.0
            };
            if distance(self.pos[p], self.pos[q]) > cutoff + slack {
                continue;
            }
            if !self.move_to(q, t) {
                continue;
            }
            self.draw_pair(p, q, t, fresh, horizon, cutoff);
        }
    }

    /// Draw the reaction time of the pair `(i, j)`, both at their positions
    /// of `t`: the earliest of its reactions and mirror images.
    fn draw_pair(&mut self, i: usize, j: usize, t: f64, fresh: bool, horizon: f64, cutoff: f64) {
        let model = self.model;
        let reactions = model.pairs_of(self.species[i] as usize, self.species[j] as usize);
        if reactions.is_empty() {
            return;
        }
        let (xi, xj) = (self.pos[i], self.pos[j]);
        // (dt, reaction, image) of the earliest reaction
        let mut best: Option<(f64, u32, [i8; 3])> = None;
        let mut look_again = f64::INFINITY;
        let real = distance(xi, xj);
        if real >= cutoff {
            return;
        }
        for &r in reactions {
            let p = &model.pairs[r as usize];
            if real < p.r_eff {
                if fresh && self.opt.start_contact && contact_reacts(&mut self.rng, p) {
                    self.push(
                        t,
                        Kind::Pair {
                            i: i as u32,
                            j: j as u32,
                            reaction: r,
                            image: [0; 3],
                            t_sample: t,
                            at: (xi, xj),
                            contact: true,
                            epoch: self.epoch,
                        },
                    );
                    return;
                }
                look_again = look_again.min(p.r_eff * p.r_eff / p.d_sum);
                continue;
            }
            self.stats.pairs_sampled += 1;
            if let Some(dt) = pair_time(&mut self.rng, p, real)
                && dt < best.map_or(horizon, |b| b.0)
            {
                best = Some((dt, r, [0; 3]));
            }
        }
        if let Boundary::Reflective { lo, hi } = self.opt.boundary
            && self.opt.wall_images
        {
            let mut choices = [[0i8; 3]; 3];
            let mut n_choices = [1usize; 3];
            for a in 0..3 {
                if (xi[a] - lo[a]) + (xj[a] - lo[a]) < cutoff {
                    choices[a][n_choices[a]] = -1;
                    n_choices[a] += 1;
                }
                if (hi[a] - xi[a]) + (hi[a] - xj[a]) < cutoff {
                    choices[a][n_choices[a]] = 1;
                    n_choices[a] += 1;
                }
            }
            for &ix in &choices[0][..n_choices[0]] {
                for &iy in &choices[1][..n_choices[1]] {
                    for &iz in &choices[2][..n_choices[2]] {
                        let image = [ix, iy, iz];
                        if image == [0; 3] {
                            continue;
                        }
                        let d = distance(xi, mirror(lo, hi, xj, image));
                        if d >= cutoff {
                            continue;
                        }
                        for &r in reactions {
                            let p = &model.pairs[r as usize];
                            if d < p.r_eff {
                                look_again = look_again.min(p.r_eff * p.r_eff / p.d_sum);
                                continue;
                            }
                            self.stats.pairs_sampled += 1;
                            if let Some(dt) = pair_time(&mut self.rng, p, d)
                                && dt < best.map_or(horizon, |b| b.0)
                            {
                                best = Some((dt, r, image));
                            }
                        }
                    }
                }
            }
        }
        match best {
            Some((dt, reaction, image)) if dt <= look_again => self.push(
                t + dt,
                Kind::Pair {
                    i: i as u32,
                    j: j as u32,
                    reaction,
                    image,
                    t_sample: t,
                    at: (xi, xj),
                    contact: false,
                    epoch: self.epoch,
                },
            ),
            _ if look_again < horizon => self.push(
                t + look_again,
                Kind::Relook {
                    i: i as u32,
                    j: j as u32,
                    epoch: self.epoch,
                },
            ),
            _ => {}
        }
    }

    /// Add the products of a reaction at `site` and draw their times.
    fn products(&mut self, products: &[usize], site: Vec3, t: f64, fresh: bool) {
        let site = match self.opt.boundary {
            Boundary::Outflow { lo, hi } if outside(lo, hi, site) => {
                for &s in products {
                    self.stats.absorbed[s] += 1;
                }
                return;
            }
            Boundary::Reflective { lo, hi } => fold(lo, hi, site),
            _ => site,
        };
        let first = self.species.len();
        for &s in products {
            self.add(s as u32, site, t);
        }
        for p in first..self.species.len() {
            self.draw(p, t, fresh);
        }
    }

    fn fire(&mut self, event: Event) {
        let model = self.model;
        let t = event.time;
        match event.kind {
            Kind::Pair {
                i,
                j,
                reaction,
                image,
                t_sample,
                at,
                contact,
                epoch,
            } => {
                let (i, j) = (i as usize, j as usize);
                if epoch != self.epoch || !self.alive[i] || !self.alive[j] {
                    return;
                }
                let p = &model.pairs[reaction as usize];
                let d = (
                    model.diffusion[self.species[i] as usize],
                    model.diffusion[self.species[j] as usize],
                );
                let partner = match self.opt.boundary {
                    Boundary::Reflective { lo, hi } if image != [0; 3] => {
                        self.stats.image_reactions += 1;
                        mirror(lo, hi, at.1, image)
                    }
                    _ => at.1,
                };
                let met =
                    reaction_positions(&mut self.rng, p.r_eff, d, (at.0, partner), t - t_sample);
                let site = product_position(&mut self.rng, d, met);
                self.kill(i);
                self.kill(j);
                self.t = t;
                if contact {
                    self.stats.contact_reactions += 1;
                } else {
                    self.stats.pair_reactions += 1;
                }
                self.stats.by_reaction[p.index] += 1;
                self.products(&p.products, site, t, contact);
            }
            Kind::Molecule { i, reaction } => {
                let i = i as usize;
                if !self.alive[i] || !self.move_to(i, t) {
                    return;
                }
                let m = &model.molecule[reaction as usize];
                let site = self.pos[i];
                self.kill(i);
                self.t = t;
                self.stats.molecule_reactions += 1;
                self.stats.by_reaction[m.index] += 1;
                if let Some(b) = m.background {
                    self.stats.background_consumed[b] += 1;
                }
                self.products(&m.products, site, t, false);
            }
            Kind::Relook { i, j, epoch } => {
                let (i, j) = (i as usize, j as usize);
                if epoch != self.epoch || !self.alive[i] || !self.alive[j] {
                    return;
                }
                self.stats.moves += 1;
                if !self.move_to(i, t) || !self.move_to(j, t) {
                    return;
                }
                self.t = t;
                let (horizon, _) = self.horizon(t);
                if horizon > 0.0 {
                    let cutoff = self.cutoff(horizon);
                    self.draw_pair(i, j, t, false, horizon, cutoff);
                }
            }
            Kind::BoxPass => {
                self.box_pass_scheduled = false;
                let (horizon, limit) = self.horizon(t);
                if horizon.is_nan() || horizon <= 0.0 {
                    return;
                }
                self.stats.box_moves += 1;
                self.stats.passes += 1;
                // Every pair time drawn so far is void.
                self.epoch += 1;
                self.sync(t);
                self.t = t;
                let cutoff = self.cutoff(horizon);
                let ids: Vec<usize> = (0..self.species.len()).filter(|i| self.alive[*i]).collect();
                let pos: Vec<Vec3> = ids.iter().map(|i| self.pos[*i]).collect();
                let mut found: Vec<(usize, usize)> = Vec::new();
                let mut grid = std::mem::take(&mut self.grid);
                grid.for_each_pair(&pos, cutoff, |a, b, _| found.push((ids[a], ids[b])));
                self.grid = grid;
                for (a, b) in found {
                    self.draw_pair(a, b, t, false, horizon, cutoff);
                }
                if let Some(limit) = limit
                    && t + limit < self.opt.t_stop
                    && ids.len() >= 2
                {
                    self.box_pass_scheduled = true;
                    self.push(t + limit, Kind::BoxPass);
                }
            }
        }
    }
}
