//! Event loop of the IRT stage with full resampling.
//!
//! Two kinds of reactions fire: reactions of a pair of particles, and
//! reactions of one particle (first order, with a fixed species such as
//! water, or with a background species). The second kind has an exponential
//! reaction time and consumes no background: the background keeps its
//! concentration, and what it loses is only counted.
//!
//! Between reflective walls every position is folded back into the box, and
//! a pair near a wall also samples a reaction time against the mirror image
//! of its partner, across every wall close enough, with the edge and corner
//! images. The image is one more independent pair, as IRT treats every
//! pair. Only the first image across each wall is used, so the look-ahead
//! of a pass is limited to the time a pair needs to feel the opposite wall;
//! past it every particle is moved and a new pass starts.
//!
//! The particles sit at their positions of the stage time `t`. One pass
//! samples a reaction time for every pair inside the cutoff radius and keeps
//! the earliest. Firing it moves the time to the reaction, puts the two
//! reactants and their products at sampled positions, moves every other
//! particle by its free diffusion over the step, and discards every sampled
//! time: the next pass draws them all again from the new positions.
//!
//! A pair found inside its effective radius reacts at once, with the contact
//! probability of the reaction, when one of its two particles has just been
//! injected (or created by such a reaction). Later, a pair that diffusion
//! brought inside its effective radius gets no reaction time in that pass:
//! its sampled time said it had not reacted yet. So that such a pair does
//! not wait for a reaction elsewhere, the next step is then at most the
//! time the pair needs to diffuse over its radius (`r_eff^2 / D`): every
//! particle is moved over it and a new pass looks again.

use crate::irt::model::IrtModel;
use crate::irt::sample::{
    Vec3, contact_reacts, diffuse, erfc_inv, pair_time, product_position, reaction_positions,
};
use crate::irt::search::CellGrid;
use crate::ssa::rng::SsaRng;

/// What happens at the faces of the box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Boundary {
    /// Free space.
    None,
    /// A particle that a move puts outside the box `[lo, hi]` is removed.
    ///
    /// This is what Geant4-DNA does at the faces of its water volume. It is
    /// not an absorbing wall: a path that leaves the box and comes back
    /// within one move is not seen, so about half as many particles are
    /// lost as at an absorbing wall when the moves are long.
    Outflow { lo: Vec3, hi: Vec3 },
    /// Reflective walls at the faces of the box `[lo, hi]`.
    Reflective { lo: Vec3, hi: Vec3 },
}

/// Settings of a stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageOptions {
    /// End of the stage [s]: no reaction is looked for beyond it.
    pub t_stop: f64,
    /// Apply the contact rule to pairs of just injected particles.
    pub start_contact: bool,
    pub boundary: Boundary,
    /// Pairs whose probability to meet before `t_stop` is below this times
    /// `r_eff / distance` are not sampled (Geant4: 0.01).
    pub cutoff_probability: f64,
    /// Fire the reactions of one particle: first order, with a fixed
    /// species, or with a background species.
    pub molecule_reactions: bool,
    /// Between reflective walls, sample the mirror images of the partner
    /// of a pair (off: a diagnostic, the walls then only fold positions).
    pub wall_images: bool,
}

/// Event counts of a stage.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StageStats {
    /// Reactions of pairs at a sampled reaction time.
    pub pair_reactions: u64,
    /// Reactions of pairs found in contact.
    pub contact_reactions: u64,
    /// Reactions of one particle: first order, with a fixed species, or
    /// with a background species.
    pub molecule_reactions: u64,
    /// Molecules of each background species consumed by these reactions.
    /// They are counted, not removed: the background keeps its
    /// concentration.
    pub background_consumed: Vec<u64>,
    /// Firings of each network reaction.
    pub by_reaction: Vec<u64>,
    /// Particles that left the box (`Boundary::Outflow`), per species.
    pub absorbed: Vec<u64>,
    /// Moves of every particle made to look again at a pair found inside
    /// its effective radius.
    pub moves: u64,
    /// Reactions of a pair with a mirror image of its partner (counted in
    /// `pair_reactions` too).
    pub image_reactions: u64,
    /// Moves of every particle made at the look-ahead limit of a box with
    /// reflective walls.
    pub box_moves: u64,
    /// Sampling passes.
    pub passes: u64,
    /// Pair reaction times drawn.
    pub pairs_sampled: u64,
}

/// The particles of a stage.
#[derive(Debug, Clone, Default)]
pub struct Particles {
    pub species: Vec<u32>,
    pub pos: Vec<Vec3>,
    /// Just injected, or created by a contact reaction.
    fresh: Vec<bool>,
    /// Particles per species.
    counts: Vec<i64>,
}

impl Particles {
    fn new(n_species: usize) -> Self {
        Self {
            counts: vec![0; n_species],
            ..Self::default()
        }
    }

    pub fn len(&self) -> usize {
        self.species.len()
    }

    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }

    fn push(&mut self, species: u32, pos: Vec3, fresh: bool) {
        self.species.push(species);
        self.pos.push(pos);
        self.fresh.push(fresh);
        self.counts[species as usize] += 1;
    }

    /// Drop the particles flagged in `dead`, keeping the order of the others.
    fn remove(&mut self, dead: &[bool]) {
        let mut keep = 0;
        for i in 0..self.len() {
            if dead[i] {
                self.counts[self.species[i] as usize] -= 1;
                continue;
            }
            self.species[keep] = self.species[i];
            self.pos[keep] = self.pos[i];
            self.fresh[keep] = self.fresh[i];
            keep += 1;
        }
        self.species.truncate(keep);
        self.pos.truncate(keep);
        self.fresh.truncate(keep);
    }
}

/// The next event of a stage.
#[derive(Debug, Clone, PartialEq)]
enum Next {
    /// No reaction before `t_stop`.
    Nothing,
    /// Reaction `reaction` (index in `IrtModel::pairs`) of the particles `i`
    /// and `j`, `dt` after the stage time.
    Pair {
        dt: f64,
        i: usize,
        j: usize,
        reaction: usize,
        /// Mirror image of `j` the reaction is with: per axis 0 (none),
        /// -1 (across the lower wall) or 1 (across the upper wall).
        image: [i8; 3],
    },
    /// No reaction within the look-ahead limit `dt` of a box with
    /// reflective walls: move every particle over `dt` and look again.
    BoxMove { dt: f64 },
    /// Reaction `reaction` (index in `IrtModel::molecule`) of the particle
    /// `i`, `dt` after the stage time.
    Molecule { dt: f64, i: usize, reaction: usize },
    /// Pairs in contact that react now: `(i, j, reaction)`.
    Contacts(Vec<(usize, usize, usize)>),
    /// No reaction within `dt`, but a pair sits inside its effective radius:
    /// move every particle over `dt` and look again.
    Move { dt: f64 },
}

/// One trajectory of the IRT stage.
#[derive(Debug, Clone)]
pub struct Stage<'a> {
    model: &'a IrtModel,
    opt: StageOptions,
    rng: SsaRng,
    particles: Particles,
    /// Time of the particle positions [s].
    t: f64,
    next: Option<Next>,
    stats: StageStats,
    grid: CellGrid,
    dead: Vec<bool>,
    /// Particles that have a reaction partner among the particles, and
    /// their positions (scratch of the pair search).
    active: Vec<u32>,
    active_pos: Vec<Vec3>,
}

/// Whether `x` is outside the box `[lo, hi]` or on one of its faces.
pub(crate) fn outside(lo: Vec3, hi: Vec3, x: Vec3) -> bool {
    (0..3).any(|a| x[a] <= lo[a] || x[a] >= hi[a])
}

/// `x` folded back into the box `[lo, hi]` by reflections at its faces.
pub fn fold(lo: Vec3, hi: Vec3, x: Vec3) -> Vec3 {
    std::array::from_fn(|a| {
        let l = hi[a] - lo[a];
        let y = (x[a] - lo[a]).rem_euclid(2.0 * l);
        lo[a] + if y > l { 2.0 * l - y } else { y }
    })
}

/// Mirror image of `x` across the walls named by `image` (per axis 0, -1 for
/// the lower wall, 1 for the upper wall).
pub(crate) fn mirror(lo: Vec3, hi: Vec3, x: Vec3, image: [i8; 3]) -> Vec3 {
    std::array::from_fn(|a| match image[a] {
        -1 => 2.0 * lo[a] - x[a],
        1 => 2.0 * hi[a] - x[a],
        _ => x[a],
    })
}

impl<'a> Stage<'a> {
    /// Empty stage at `t_start`.
    pub fn new(model: &'a IrtModel, opt: StageOptions, rng: SsaRng, t_start: f64) -> Self {
        let n = model.n_species();
        Self {
            model,
            opt,
            rng,
            particles: Particles::new(n),
            t: t_start,
            next: None,
            stats: StageStats {
                by_reaction: vec![0; model.n_reactions()],
                absorbed: vec![0; n],
                background_consumed: vec![0; n],
                ..StageStats::default()
            },
            grid: CellGrid::new(),
            dead: Vec::new(),
            active: Vec::new(),
            active_pos: Vec::new(),
        }
    }

    /// Time of the particle positions [s].
    pub fn time(&self) -> f64 {
        self.t
    }

    /// Particles per species.
    pub fn counts(&self) -> &[i64] {
        &self.particles.counts
    }

    pub fn particles(&self) -> &Particles {
        &self.particles
    }

    /// `(species, position)` of every particle, at the stage time.
    pub fn molecules(&self) -> Vec<(u32, Vec3)> {
        self.particles
            .species
            .iter()
            .copied()
            .zip(self.particles.pos.iter().copied())
            .collect()
    }

    pub fn stats(&self) -> &StageStats {
        &self.stats
    }

    /// Add `(species, position)` particles at the stage time. Call
    /// [`Stage::sync`] first when the stage is behind the injection time.
    pub fn inject(&mut self, molecules: &[(u32, Vec3)]) {
        for &(s, pos) in molecules {
            self.particles.push(s, pos, self.opt.start_contact);
        }
        self.next = None;
    }

    /// Fire every reaction up to `t`. The particles stay at the positions of
    /// the last reaction: only their counts are those of `t`.
    pub fn advance_to(&mut self, t: f64) {
        loop {
            if self.next.is_none() {
                self.next = Some(self.sample());
            }
            match self.next.take() {
                Some(Next::Contacts(list)) => self.fire_contacts(&list),
                Some(Next::Pair {
                    dt,
                    i,
                    j,
                    reaction,
                    image,
                }) if self.t + dt <= t => self.fire_pair(dt, i, j, reaction, image),
                Some(Next::BoxMove { dt }) if self.t + dt <= t => {
                    self.stats.box_moves += 1;
                    self.sync(self.t + dt);
                }
                Some(Next::Molecule { dt, i, reaction }) if self.t + dt <= t => {
                    self.fire_molecule(dt, i, reaction)
                }
                Some(Next::Move { dt }) if self.t + dt <= t => {
                    self.stats.moves += 1;
                    self.sync(self.t + dt);
                }
                other => {
                    self.next = other;
                    return;
                }
            }
        }
    }

    /// Move every particle to `t` by its free diffusion and forget the
    /// sampled times. Reactions up to `t` must have fired
    /// ([`Stage::advance_to`]).
    pub fn sync(&mut self, t: f64) {
        let dt = t - self.t;
        if dt <= 0.0 {
            return;
        }
        self.dead.clear();
        self.dead.resize(self.particles.len(), false);
        self.diffuse_all(dt, usize::MAX, usize::MAX);
        self.particles.remove(&self.dead);
        self.t = t;
        self.next = None;
    }

    /// Cutoff radius of the pair search for the time left to `t_stop`.
    fn cutoff(&self, horizon: f64) -> f64 {
        self.model.max_radius()
            + 2.0
                * (self.model.max_d_sum() * horizon).sqrt()
                * erfc_inv(self.opt.cutoff_probability)
    }

    /// One sampling pass: the earliest reaction, or the contact reactions.
    fn sample(&mut self) -> Next {
        let horizon = self.opt.t_stop - self.t;
        if horizon.is_nan() || horizon <= 0.0 || self.particles.is_empty() {
            return Next::Nothing;
        }
        self.stats.passes += 1;
        // Between reflective walls a pass looks ahead no further than the
        // time for which the first mirror images are enough: its cutoff
        // radius is then the shortest edge of the box.
        let full_horizon = horizon;
        let horizon = match self.opt.boundary {
            Boundary::Reflective { lo, hi } if self.model.max_d_sum() > 0.0 => {
                let edge = (0..3).map(|a| hi[a] - lo[a]).fold(f64::INFINITY, f64::min);
                let reach = edge / (2.0 * erfc_inv(self.opt.cutoff_probability));
                horizon.min(reach * reach / self.model.max_d_sum())
            }
            _ => horizon,
        };
        let cutoff = self.cutoff(horizon);
        let model = self.model;
        let Self {
            particles,
            rng,
            grid,
            stats,
            opt,
            active,
            active_pos,
            ..
        } = self;
        let species = &particles.species;
        let fresh = &particles.fresh;
        // Only the particles of a species that has a partner present can
        // react in a pair: the others are left out of the search.
        let mut partner = vec![false; model.n_species()];
        for p in &model.pairs {
            let present = if p.a == p.b {
                particles.counts[p.a] >= 2
            } else {
                particles.counts[p.a] >= 1 && particles.counts[p.b] >= 1
            };
            if present {
                partner[p.a] = true;
                partner[p.b] = true;
            }
        }
        active.clear();
        active_pos.clear();
        for (i, s) in species.iter().enumerate() {
            if partner[*s as usize] {
                active.push(i as u32);
                active_pos.push(particles.pos[i]);
            }
        }
        let mut best = Next::Nothing;
        let mut best_dt = horizon;
        let mut contacts: Vec<(usize, usize, usize)> = Vec::new();
        // Time for the closest overlapping pair to diffuse over its radius.
        let mut look_again = f64::INFINITY;
        grid.for_each_pair(active_pos, cutoff, |i, j, d2| {
            let (i, j) = (active[i] as usize, active[j] as usize);
            let reactions = model.pairs_of(species[i] as usize, species[j] as usize);
            if reactions.is_empty() {
                return;
            }
            let distance = d2.sqrt();
            for &r in reactions {
                let p = &model.pairs[r as usize];
                if d2 < p.r_eff * p.r_eff {
                    if opt.start_contact && (fresh[i] || fresh[j]) && contact_reacts(rng, p) {
                        contacts.push((i, j, r as usize));
                        break;
                    }
                    look_again = look_again.min(p.r_eff * p.r_eff / p.d_sum);
                    continue;
                }
                stats.pairs_sampled += 1;
                if let Some(dt) = pair_time(rng, p, distance)
                    && dt < best_dt
                {
                    best_dt = dt;
                    best = Next::Pair {
                        dt,
                        i,
                        j,
                        reaction: r as usize,
                        image: [0; 3],
                    };
                }
            }
            // The mirror images of `j` across the walls both particles are
            // close to: an image is further than `j` itself, so every image
            // pair inside the cutoff is among the pairs found here.
            let Boundary::Reflective { lo, hi } = opt.boundary else {
                return;
            };
            if !opt.wall_images {
                return;
            }
            let (xi, xj) = (particles.pos[i], particles.pos[j]);
            // Per axis: the images worth looking at (0 is always one).
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
                        let m = mirror(lo, hi, xj, image);
                        let d2 = (0..3).map(|a| (xi[a] - m[a]).powi(2)).sum::<f64>();
                        if d2 >= cutoff * cutoff {
                            continue;
                        }
                        let distance = d2.sqrt();
                        for &r in reactions {
                            let p = &model.pairs[r as usize];
                            if d2 < p.r_eff * p.r_eff {
                                look_again = look_again.min(p.r_eff * p.r_eff / p.d_sum);
                                continue;
                            }
                            stats.pairs_sampled += 1;
                            if let Some(dt) = pair_time(rng, p, distance)
                                && dt < best_dt
                            {
                                best_dt = dt;
                                best = Next::Pair {
                                    dt,
                                    i,
                                    j,
                                    reaction: r as usize,
                                    image,
                                };
                            }
                        }
                    }
                }
            }
        });
        let n_active = active.len();
        particles.fresh.fill(false);
        if !contacts.is_empty() {
            return Next::Contacts(contacts);
        }
        // The reactions of single particles: one exponential time for the
        // sum of their rates, then the particle and its reaction.
        if opt.molecule_reactions {
            let weight = |s: usize| particles.counts[s] as f64 * model.molecule_rate(s);
            let total: f64 = (0..model.n_species()).map(weight).sum();
            if total > 0.0 {
                let dt = rng.exp1() / total;
                if dt < best_dt && dt < look_again {
                    let mut u = rng.uniform() * total;
                    let mut chosen = 0;
                    for s in 0..model.n_species() {
                        let w = weight(s);
                        if w > 0.0 {
                            chosen = s;
                            if u < w {
                                break;
                            }
                            u -= w;
                        }
                    }
                    let k = rng.below(particles.counts[chosen] as u64) as usize;
                    let i = (0..species.len())
                        .filter(|i| species[*i] as usize == chosen)
                        .nth(k)
                        .expect("a particle of the chosen species");
                    let list = model.molecule_of(chosen);
                    let mut v = rng.uniform() * model.molecule_rate(chosen);
                    let mut reaction = list[list.len() - 1] as usize;
                    for &r in list {
                        let rate = model.molecule[r as usize].rate;
                        if v < rate {
                            reaction = r as usize;
                            break;
                        }
                        v -= rate;
                    }
                    return Next::Molecule { dt, i, reaction };
                }
            }
        }
        if look_again < best_dt {
            Next::Move { dt: look_again }
        } else if best == Next::Nothing && horizon < full_horizon && n_active >= 2 {
            Next::BoxMove { dt: horizon }
        } else {
            best
        }
    }

    /// Move every particle but `skip_a` and `skip_b` over `dt`, flagging in
    /// `dead` those that leave the box.
    fn diffuse_all(&mut self, dt: f64, skip_a: usize, skip_b: usize) {
        for k in 0..self.particles.len() {
            if k == skip_a || k == skip_b {
                continue;
            }
            let s = self.particles.species[k] as usize;
            let d = self.model.diffusion[s];
            diffuse(&mut self.rng, d, dt, &mut self.particles.pos[k]);
            match self.opt.boundary {
                Boundary::Outflow { lo, hi } if outside(lo, hi, self.particles.pos[k]) => {
                    self.dead[k] = true;
                    self.stats.absorbed[s] += 1;
                }
                Boundary::Reflective { lo, hi } => {
                    self.particles.pos[k] = fold(lo, hi, self.particles.pos[k]);
                }
                _ => {}
            }
        }
    }

    /// Add the products of `reaction` at `site`.
    fn add_products(&mut self, reaction: usize, site: Vec3, fresh: bool) {
        let model = self.model;
        if let Boundary::Outflow { lo, hi } = self.opt.boundary
            && outside(lo, hi, site)
        {
            for &s in &model.pairs[reaction].products {
                self.stats.absorbed[s] += 1;
            }
            return;
        }
        let site = match self.opt.boundary {
            Boundary::Reflective { lo, hi } => fold(lo, hi, site),
            _ => site,
        };
        for &s in &model.pairs[reaction].products {
            self.particles.push(s as u32, site, fresh);
        }
    }

    fn fire_pair(&mut self, dt: f64, i: usize, j: usize, reaction: usize, image: [i8; 3]) {
        let model = self.model;
        let p = &model.pairs[reaction];
        let d = (
            model.diffusion[self.particles.species[i] as usize],
            model.diffusion[self.particles.species[j] as usize],
        );
        // A reaction with a mirror image happens where the image is: the
        // products are folded back into the box.
        let partner = match self.opt.boundary {
            Boundary::Reflective { lo, hi } if image != [0; 3] => {
                self.stats.image_reactions += 1;
                mirror(lo, hi, self.particles.pos[j], image)
            }
            _ => self.particles.pos[j],
        };
        let at = (self.particles.pos[i], partner);
        let met = reaction_positions(&mut self.rng, p.r_eff, d, at, dt);
        let site = product_position(&mut self.rng, d, met);
        self.dead.clear();
        self.dead.resize(self.particles.len(), false);
        self.diffuse_all(dt, i, j);
        self.dead[i] = true;
        self.dead[j] = true;
        self.particles.remove(&self.dead);
        self.t += dt;
        self.add_products(reaction, site, false);
        self.stats.pair_reactions += 1;
        self.stats.by_reaction[p.index] += 1;
    }

    /// React the particle `i` alone, `dt` after the stage time: every
    /// particle moves over `dt`, then the products replace `i` where it is.
    fn fire_molecule(&mut self, dt: f64, i: usize, reaction: usize) {
        let model = self.model;
        let m = &model.molecule[reaction];
        self.dead.clear();
        self.dead.resize(self.particles.len(), false);
        self.diffuse_all(dt, usize::MAX, usize::MAX);
        self.t += dt;
        // A particle that left the box on the way does not react.
        if self.dead[i] {
            self.particles.remove(&self.dead);
            return;
        }
        let site = self.particles.pos[i];
        self.dead[i] = true;
        self.particles.remove(&self.dead);
        for &s in &m.products {
            self.particles.push(s as u32, site, false);
        }
        self.stats.molecule_reactions += 1;
        self.stats.by_reaction[m.index] += 1;
        if let Some(b) = m.background {
            self.stats.background_consumed[b] += 1;
        }
    }

    /// React the pairs in contact, in the order they were found; a pair one
    /// of whose particles has already reacted is skipped.
    fn fire_contacts(&mut self, list: &[(usize, usize, usize)]) {
        let model = self.model;
        self.dead.clear();
        self.dead.resize(self.particles.len(), false);
        for &(i, j, reaction) in list {
            if self.dead[i] || self.dead[j] {
                continue;
            }
            let d = (
                model.diffusion[self.particles.species[i] as usize],
                model.diffusion[self.particles.species[j] as usize],
            );
            let at = (self.particles.pos[i], self.particles.pos[j]);
            let site = product_position(&mut self.rng, d, at);
            self.dead[i] = true;
            self.dead[j] = true;
            self.add_products(reaction, site, true);
            self.stats.contact_reactions += 1;
            self.stats.by_reaction[model.pairs[reaction].index] += 1;
        }
        self.dead.resize(self.particles.len(), false);
        self.particles.remove(&self.dead);
    }
}
