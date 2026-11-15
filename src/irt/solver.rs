//! Event loop of the IRT stage with full resampling.
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
}

/// Event counts of a stage.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StageStats {
    /// Reactions of pairs at a sampled reaction time.
    pub pair_reactions: u64,
    /// Reactions of pairs found in contact.
    pub contact_reactions: u64,
    /// Firings of each network reaction.
    pub by_reaction: Vec<u64>,
    /// Particles that left the box (`Boundary::Outflow`), per species.
    pub absorbed: Vec<u64>,
    /// Moves of every particle made to look again at a pair found inside
    /// its effective radius.
    pub moves: u64,
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
    },
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
}

/// Whether `x` is outside the box `[lo, hi]` or on one of its faces.
fn outside(lo: Vec3, hi: Vec3, x: Vec3) -> bool {
    (0..3).any(|a| x[a] <= lo[a] || x[a] >= hi[a])
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
                ..StageStats::default()
            },
            grid: CellGrid::new(),
            dead: Vec::new(),
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
                Some(Next::Pair { dt, i, j, reaction }) if self.t + dt <= t => {
                    self.fire_pair(dt, i, j, reaction)
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
        if horizon.is_nan() || horizon <= 0.0 || self.particles.len() < 2 {
            return Next::Nothing;
        }
        self.stats.passes += 1;
        let cutoff = self.cutoff(horizon);
        let model = self.model;
        let Self {
            particles,
            rng,
            grid,
            stats,
            opt,
            ..
        } = self;
        let species = &particles.species;
        let fresh = &particles.fresh;
        let mut best = Next::Nothing;
        let mut best_dt = horizon;
        let mut contacts: Vec<(usize, usize, usize)> = Vec::new();
        // Time for the closest overlapping pair to diffuse over its radius.
        let mut look_again = f64::INFINITY;
        grid.for_each_pair(&particles.pos, cutoff, |i, j, d2| {
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
                    };
                }
            }
        });
        particles.fresh.fill(false);
        if !contacts.is_empty() {
            Next::Contacts(contacts)
        } else if look_again < best_dt {
            Next::Move { dt: look_again }
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
            if let Boundary::Outflow { lo, hi } = self.opt.boundary
                && outside(lo, hi, self.particles.pos[k])
            {
                self.dead[k] = true;
                self.stats.absorbed[s] += 1;
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
        for &s in &model.pairs[reaction].products {
            self.particles.push(s as u32, site, fresh);
        }
    }

    fn fire_pair(&mut self, dt: f64, i: usize, j: usize, reaction: usize) {
        let model = self.model;
        let p = &model.pairs[reaction];
        let d = (
            model.diffusion[self.particles.species[i] as usize],
            model.diffusion[self.particles.species[j] as usize],
        );
        let at = (self.particles.pos[i], self.particles.pos[j]);
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
