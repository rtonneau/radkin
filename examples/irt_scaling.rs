//! Run time of the IRT stage against the number of overlapping tracks, for
//! the two pair rules (`docs/perf/irt-scaling.md`).
//!
//! The tracks are the ten reference events of `tests/data/irt/` (2 um long
//! segments of 100 keV electron tracks, 87 molecules each on average), used
//! in turn, injected together at 1.001 ps and moved sideways by a random
//! offset inside a square of `--spread-nm`: they overlap. The chemistry is
//! `data/chem_g4.toml` with its background reactions, in free space. One
//! trajectory is timed per case.
//!
//! ```text
//! cargo run --release --example irt_scaling -- --tracks 1,10,100
//! cargo run --release --example irt_scaling -- --tracks 100,300 --full-up-to 0
//! ```

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use radkin::irt::{IrtSimulation, run_irt_with};
use radkin::ssa::SsaRng;

#[derive(Parser)]
struct Args {
    /// Numbers of tracks, comma separated.
    #[arg(long, value_delimiter = ',', default_value = "1,3,10,30,100")]
    tracks: Vec<usize>,
    /// End of the IRT stage [s].
    #[arg(long, default_value_t = 1e-6)]
    t_end: f64,
    /// Edge of the square the track axes are spread over [nm].
    #[arg(long, default_value_t = 1000.0)]
    spread_nm: f64,
    /// Largest number of tracks run with full resampling (its cost grows with
    /// the square of the molecule count per reaction).
    #[arg(long, default_value_t = 10)]
    full_up_to: usize,
    /// Random seed of the offsets and of the trajectory.
    #[arg(long, default_value_t = 1)]
    seed: u64,
}

fn toml_path(p: &std::path::Path) -> String {
    p.display().to_string().replace('\\', "/")
}

fn main() -> radkin::Result<()> {
    let args = Args::parse();
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let molecules = toml_path(&root.join("tests/data/irt/g4_events_molecules.csv"));
    let dir = std::env::temp_dir().join(format!("radkin_irt_scaling_{}", std::process::id()));
    std::fs::create_dir_all(&dir).map_err(|e| radkin::RadkinError::io(&dir, e))?;
    println!("tracks,molecules,resampling,wall_s,reactions,pair_times_drawn,e_aq_left");
    for &n in &args.tracks {
        let mut rng = SsaRng::new(args.seed, 1000 + n as u64);
        let mut seeds = String::new();
        for k in 0..n {
            let x = args.spread_nm * (rng.uniform() - 0.5);
            let y = args.spread_nm * (rng.uniform() - 0.5);
            seeds.push_str(&format!(
                "[[irt.seed]]\nfile = \"{molecules}\"\ntime = 1.001e-12\nevents = [{}]\noffset_nm = [{x}, {y}, 0]\n",
                k % 10
            ));
        }
        for rule in ["incremental", "full"] {
            if rule == "full" && n > args.full_up_to {
                continue;
            }
            let text = format!(
                "[includes]\ntoml = [\"radkin/chem_g4.toml\"]\n[irt]\nresampling = \"{rule}\"\n\
                 start_contact = false\nt_end = {:e}\nrecord_times = []\n{seeds}",
                args.t_end
            );
            let path = dir.join(format!("{rule}_{n}.toml"));
            std::fs::write(&path, text).map_err(|e| radkin::RadkinError::io(&path, e))?;
            let sim = IrtSimulation::from_toml(&path)?;
            let n_molecules: usize = sim.seeds.iter().map(|s| s.molecules.len()).sum();
            let started = Instant::now();
            let out = run_irt_with(&sim, 1, args.seed)?;
            let wall = started.elapsed().as_secs_f64();
            let e_aq = out.species.iter().position(|s| s == "e_aq").unwrap_or(0);
            println!(
                "{n},{n_molecules},{rule},{wall:.3},{},{},{}",
                out.stats.pair_reactions
                    + out.stats.contact_reactions
                    + out.stats.molecule_reactions,
                out.stats.pairs_sampled,
                out.totals_mean.last().map_or(0.0, |row| row[e_aq])
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    Ok(())
}
