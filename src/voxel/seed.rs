//! Seeds: Geant4-DNA molecules binned into voxel concentrations and injected
//! at given times.
//!
//! Molecule CSV (header required, `#` comments allowed):
//!
//! ```text
//! species,x,y,z,event_id
//! e_aq^-1,1.5,-2.0,30.1,0
//! °OH^0,1.7,-2.2,30.0,0
//! ```
//!
//! Positions are in nm; `event_id` is optional. Species names go through the
//! alias table ([`DEFAULT_GEANT4_ALIASES`] plus `[voxel.geant4.aliases]`); a
//! name that is not an alias is used as a radkin label as is.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;

use crate::chemistry::species::ListSpecies;
use crate::error::{RadkinError, Result};
use crate::voxel::driver::VoxelRun;
use crate::voxel::grid::{VoxelGrid, VoxelState};

/// Geant4-DNA 11.4.3 molecule names -> radkin labels (option3 /
/// `G4ChemDissociationChannels_option1::ConstructMolecule`). Both the
/// `G4MolecularConfiguration::GetName()` form (`<definition>^<charge>`) and
/// the configuration tags created there are listed.
pub const DEFAULT_GEANT4_ALIASES: [(&str, &str); 30] = [
    ("e_aq^-1", "e_aq"),
    ("°OH^0", "OH_r"),
    ("H3O^1", "H3O_plus"),
    ("H^0", "H_r"),
    ("H_2^0", "H2"),
    ("H2O2^0", "H2O2"),
    ("OH^-1", "OH_minus"),
    ("HO_2°^0", "HO2"),
    ("HO_2^-1", "HO2_minus"),
    ("O_2^0", "O2"),
    ("O_2^-1", "O2_minus"),
    ("O_3^0", "O3"),
    ("O_3^-1", "O3_minus"),
    ("O^-1", "O_minus"),
    ("°O^0", "O_atom"),
    ("e_aq", "e_aq"),
    ("°OH", "OH_r"),
    ("H3Op", "H3O_plus"),
    ("H", "H_r"),
    ("H2", "H2"),
    ("H2O2", "H2O2"),
    ("OHm", "OH_minus"),
    ("HO2°", "HO2"),
    ("HO2m", "HO2_minus"),
    ("O2", "O2"),
    ("O2m", "O2_minus"),
    ("O3", "O3"),
    ("O3m", "O3_minus"),
    ("Om", "O_minus"),
    ("Oxy", "O_atom"),
];

/// One molecule of a seed file.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    pub species: String,
    /// Position [m].
    pub pos: [f64; 3],
    pub event_id: Option<i64>,
}

/// Read a molecule CSV (positions converted from nm to m).
pub fn read_molecules(path: &Path) -> Result<Vec<Molecule>> {
    let text = std::fs::read_to_string(path).map_err(|e| RadkinError::io(path, e))?;
    parse_molecules(&text).map_err(|e| RadkinError::Config(format!("{}: {e}", path.display())))
}

pub fn parse_molecules(text: &str) -> Result<Vec<Molecule>> {
    let mut lines = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.trim()))
        .filter(|(_, l)| !l.is_empty() && !l.starts_with('#'));
    let (_, header) = lines
        .next()
        .ok_or_else(|| RadkinError::Config("empty molecule file".into()))?;
    let cols: Vec<String> = header
        .split(',')
        .map(|c| c.trim().to_ascii_lowercase())
        .collect();
    let col = |name: &str| cols.iter().position(|c| c == name);
    let need = |name: &str| {
        col(name).ok_or_else(|| RadkinError::Config(format!("missing column '{name}'")))
    };
    let (is, ix, iy, iz) = (need("species")?, need("x")?, need("y")?, need("z")?);
    let ie = col("event_id");
    let mut out = Vec::new();
    for (n, line) in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() != cols.len() {
            return Err(RadkinError::Config(format!(
                "line {n}: {} fields, expected {}",
                f.len(),
                cols.len()
            )));
        }
        let num = |i: usize| -> Result<f64> {
            f[i].parse::<f64>()
                .map_err(|_| RadkinError::Config(format!("line {n}: bad number '{}'", f[i])))
        };
        let event_id =
            match ie {
                None => None,
                Some(i) => Some(f[i].parse::<i64>().map_err(|_| {
                    RadkinError::Config(format!("line {n}: bad event_id '{}'", f[i]))
                })?),
            };
        out.push(Molecule {
            species: f[is].to_string(),
            pos: [num(ix)? * 1e-9, num(iy)? * 1e-9, num(iz)? * 1e-9],
            event_id,
        });
    }
    Ok(out)
}

/// Name mapping and out-of-grid policy shared by all seeds.
#[derive(Debug, Clone, PartialEq)]
pub struct SeedPolicy {
    pub aliases: IndexMap<String, String>,
    /// Geant4 names dropped silently (e.g. `H2O^0`).
    pub ignore: Vec<String>,
    /// Drop molecules outside the grid instead of failing.
    pub clip: bool,
}

impl Default for SeedPolicy {
    fn default() -> Self {
        Self {
            aliases: DEFAULT_GEANT4_ALIASES
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect(),
            ignore: Vec::new(),
            clip: false,
        }
    }
}

/// `[[voxel.seed]]` entry.
#[derive(Debug, Clone, PartialEq)]
pub struct SeedSpec {
    pub file: PathBuf,
    /// Injection time [s].
    pub time: f64,
    /// Keep only these event ids (all when `None`).
    pub events: Option<Vec<i64>>,
    /// Added to every position [m].
    pub offset: [f64; 3],
    /// Energy deposited by the seed's track(s) [eV], for G-values.
    pub energy_ev: Option<f64>,
}

/// What binning a seed did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SeedReport {
    pub placed: usize,
    pub clipped: usize,
    pub ignored: usize,
    /// Molecules placed per radkin label.
    pub per_species: BTreeMap<String, usize>,
}

/// Bin `molecules` (offset added) into an amount state [µmol/L] of the
/// dynamic species.
pub fn bin_molecules(
    molecules: &[Molecule],
    offset: [f64; 3],
    grid: &VoxelGrid,
    species: &ListSpecies,
    policy: &SeedPolicy,
) -> Result<(VoxelState, SeedReport)> {
    let mut st = VoxelState::zeros(species.n_dyn(), grid.n_voxels());
    let mut rep = SeedReport::default();
    let c = grid.molecule_concentration() * 1e6;
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for m in molecules {
        if policy.ignore.contains(&m.species) {
            rep.ignored += 1;
            continue;
        }
        let label = policy
            .aliases
            .get(&m.species)
            .map(String::as_str)
            .unwrap_or(&m.species);
        if !species.is_dyn(label) {
            return Err(RadkinError::Config(format!(
                "seed molecule '{}' (radkin '{label}') is not a dynamic species of the chemistry; add an alias or list it in voxel.geant4.ignore",
                m.species
            )));
        }
        let pos: [f64; 3] = std::array::from_fn(|a| m.pos[a] + offset[a]);
        let Some(v) = grid.locate(pos) else {
            rep.clipped += 1;
            for a in 0..3 {
                lo[a] = lo[a].min(pos[a]);
                hi[a] = hi[a].max(pos[a]);
            }
            continue;
        };
        let sp = species.index_of(label)?;
        st.field_mut(sp)[v] += c;
        rep.placed += 1;
        *rep.per_species.entry(label.to_string()).or_default() += 1;
    }
    if rep.clipped > 0 {
        let nm = |p: [f64; 3]| p.map(|x| x * 1e9);
        let (glo, ghi) = grid.bounds();
        let msg = format!(
            "{} of {} seed molecules outside the grid [{:?}, {:?}] nm (their bounding box: [{:?}, {:?}] nm)",
            rep.clipped,
            rep.clipped + rep.placed,
            nm(glo),
            nm(ghi),
            nm(lo),
            nm(hi)
        );
        if !policy.clip {
            return Err(RadkinError::Config(format!(
                "{msg}; enlarge the grid or set voxel.clip = true"
            )));
        }
        log::warn!(
            "{msg}: dropped ({:.3}% lost)",
            100.0 * rep.clipped as f64 / (rep.clipped + rep.placed) as f64
        );
    }
    Ok((st, rep))
}

/// A seed ready to inject.
#[derive(Debug, Clone)]
pub struct LoadedSeed {
    pub spec: SeedSpec,
    pub amount: VoxelState,
    pub report: SeedReport,
}

/// Read and bin every seed.
pub fn load_seeds(
    specs: &[SeedSpec],
    grid: &VoxelGrid,
    species: &ListSpecies,
    policy: &SeedPolicy,
) -> Result<Vec<LoadedSeed>> {
    specs
        .iter()
        .map(|spec| {
            let mut mols = read_molecules(&spec.file)?;
            if let Some(ev) = &spec.events {
                mols.retain(|m| m.event_id.is_some_and(|e| ev.contains(&e)));
            }
            let (amount, report) = bin_molecules(&mols, spec.offset, grid, species, policy)?;
            Ok(LoadedSeed {
                spec: spec.clone(),
                amount,
                report,
            })
        })
        .collect()
}

/// Run to `t_end`, injecting each seed at its time and calling `on_output`
/// at every `output_times` entry (and after each injection at that time).
pub fn run_schedule(
    run: &mut VoxelRun,
    seeds: &[LoadedSeed],
    output_times: &[f64],
    t_end: f64,
    mut on_output: impl FnMut(&VoxelRun) -> Result<()>,
) -> Result<()> {
    let mut seeds: Vec<&LoadedSeed> = seeds.iter().collect();
    seeds.sort_by(|a, b| a.spec.time.total_cmp(&b.spec.time));
    if let Some(s) = seeds
        .iter()
        .find(|s| s.spec.time < run.t || s.spec.time > t_end)
    {
        return Err(RadkinError::Config(format!(
            "seed time {} s is outside the run [{}, {t_end}] s",
            s.spec.time, run.t
        )));
    }
    let mut outputs: Vec<f64> = output_times
        .iter()
        .copied()
        .filter(|t| *t >= run.t && *t <= t_end)
        .collect();
    outputs.sort_by(f64::total_cmp);
    outputs.dedup();
    let (mut is, mut io) = (0, 0);
    loop {
        let next_seed = seeds.get(is).map(|s| s.spec.time);
        let next_out = outputs.get(io).copied();
        let next = match (next_seed, next_out) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => break,
        };
        run.advance_to(next)?;
        while is < seeds.len() && seeds[is].spec.time == next {
            run.inject(&seeds[is].amount);
            is += 1;
        }
        if io < outputs.len() && outputs[io] == next {
            on_output(run)?;
            io += 1;
        }
    }
    run.advance_to(t_end)
}

#[cfg(test)]
mod tests {
    /// Every option3 configuration (GetName and tag forms) is either aliased
    /// to a radkin label or a documented non-solute.
    #[test]
    fn option3_names_all_covered() {
        // (GetName, tag) of G4ChemDissociationChannels_option1::ConstructMolecule.
        const SOLUTES: [(&str, &str); 15] = [
            ("H3O^1", "H3Op"),
            ("°OH^0", "°OH"),
            ("OH^-1", "OHm"),
            ("e_aq^-1", "e_aq"),
            ("H^0", "H"),
            ("H_2^0", "H2"),
            ("H2O2^0", "H2O2"),
            ("HO_2°^0", "HO2°"),
            ("HO_2^-1", "HO2m"),
            ("°O^0", "Oxy"),
            ("O^-1", "Om"),
            ("O_2^0", "O2"),
            ("O_2^-1", "O2m"),
            ("O_3^0", "O3"),
            ("O_3^-1", "O3m"),
        ];
        // Water and bulk configurations, not solutes: no alias on purpose.
        const NON_SOLUTE: [&str; 5] = ["H2O^0", "H2O(B)", "H3Op(B)", "OHm(B)", "NoneM"];
        let table: std::collections::HashMap<_, _> =
            DEFAULT_GEANT4_ALIASES.iter().copied().collect();
        for (name, tag) in SOLUTES {
            let a = table
                .get(name)
                .unwrap_or_else(|| panic!("no alias for {name}"));
            let b = table
                .get(tag)
                .unwrap_or_else(|| panic!("no alias for {tag}"));
            assert_eq!(a, b, "{name} vs {tag}");
        }
        for n in NON_SOLUTE {
            assert!(!table.contains_key(n), "{n} must not be aliased");
        }
        assert_eq!(SOLUTES.len() + NON_SOLUTE.len(), 20);
    }

    use super::*;
    use crate::chemistry::env::SimEnv;
    use crate::matrices::SimMatrices;
    use crate::voxel::diffusion::Diffusion;
    use crate::voxel::driver::StepControl;
    use crate::voxel::reaction::tests::{opts, water_env};

    const CSV: &str = "# test track\nspecies,x,y,z,event_id\ne_aq^-1,0.5,0.5,0.5,2\n°OH^0,1.5,0.5,0.5,0\nH3Op,1.6,0.4,0.5,0\nH2O^0,0,0,0,0\nH2,8.5,0.5,0.5,1\n";

    fn grid() -> VoxelGrid {
        VoxelGrid::new([10, 2, 2], [1e-9; 3], [0.0; 3]).unwrap()
    }

    fn policy() -> SeedPolicy {
        SeedPolicy {
            ignore: vec!["H2O^0".into()],
            ..SeedPolicy::default()
        }
    }

    #[test]
    fn parses_and_bins_preserving_counts() {
        let env = water_env();
        let mols = parse_molecules(CSV).unwrap();
        assert_eq!(mols.len(), 5);
        assert_eq!(mols[1].species, "°OH^0");
        let g = grid();
        let (st, rep) = bin_molecules(&mols, [0.0; 3], &g, &env.species, &policy()).unwrap();
        assert_eq!((rep.placed, rep.ignored, rep.clipped), (4, 1, 0));
        let c = g.molecule_concentration() * 1e6;
        let total: f64 = st.data.iter().sum();
        assert!((total / c - 4.0).abs() < 1e-12);
        let ioh = env.species.index_of("OH_r").unwrap();
        assert!((st.field(ioh)[g.index(1, 0, 0)] / c - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unknown_species_is_an_error() {
        let env = water_env();
        let mols = parse_molecules(CSV).unwrap();
        let err = bin_molecules(
            &mols,
            [0.0; 3],
            &grid(),
            &env.species,
            &SeedPolicy::default(),
        )
        .unwrap_err();
        assert!(err.to_string().contains("H2O^0"), "{err}");
    }

    #[test]
    fn out_of_grid_error_or_clip() {
        let env = water_env();
        let mols = parse_molecules(CSV).unwrap();
        let off = [5e-9, 0.0, 0.0];
        let err = bin_molecules(&mols, off, &grid(), &env.species, &policy()).unwrap_err();
        assert!(err.to_string().contains("1 of 4"), "{err}");
        let p = SeedPolicy {
            clip: true,
            ..policy()
        };
        let (_, rep) = bin_molecules(&mols, off, &grid(), &env.species, &p).unwrap();
        assert_eq!((rep.placed, rep.clipped), (3, 1));
    }

    #[test]
    fn bad_csv_is_an_error() {
        assert!(parse_molecules("species,x,y\nA,1,2\n").is_err());
        assert!(parse_molecules("species,x,y,z\nA,1,2\n").is_err());
        assert!(parse_molecules("species,x,y,z\nA,1,b,3\n").is_err());
    }

    #[test]
    fn two_seeds_injected_at_their_times() {
        let env: SimEnv = water_env();
        let mat = SimMatrices::build(&env).unwrap();
        let g = grid();
        let n = env.species.n_dyn();
        let ih2 = env.species.index_of("H2").unwrap();
        let mols = parse_molecules(CSV).unwrap();
        let mk = |time: f64, event: i64| {
            let m: Vec<Molecule> = mols
                .iter()
                .filter(|m| m.event_id == Some(event))
                .cloned()
                .collect();
            let (amount, report) =
                bin_molecules(&m, [0.0; 3], &g, &env.species, &policy()).unwrap();
            LoadedSeed {
                spec: SeedSpec {
                    file: PathBuf::new(),
                    time,
                    events: Some(vec![event]),
                    offset: [0.0; 3],
                    energy_ev: None,
                },
                amount,
                report,
            }
        };
        let seeds = vec![mk(1e-9, 1), mk(1e-12, 0)];
        let mut coeffs = vec![0.0; n];
        coeffs[ih2] = 5e-9;
        let mut run = VoxelRun::new(
            mat,
            Diffusion::new(g.clone(), coeffs),
            VoxelState::zeros(n, g.n_voxels()),
            0.0,
            StepControl::default(),
            opts(),
        );
        let c = g.molecule_concentration() * 1e6;
        let mut seen = Vec::new();
        run_schedule(&mut run, &seeds, &[5e-10, 1e-9, 2e-9], 2e-9, |r| {
            seen.push((r.t, r.state.field(ih2).iter().sum::<f64>() / c));
            Ok(())
        })
        .unwrap();
        assert_eq!(run.t, 2e-9);
        assert_eq!(seen.len(), 3);
        assert_eq!(seen[0].0, 5e-10);
        assert!(seen[0].1.abs() < 1e-9, "H2 before its seed: {}", seen[0].1);
        assert!(
            (seen[1].1 - 1.0).abs() < 1e-3,
            "H2 at injection: {}",
            seen[1].1
        );
    }
}
