//! The sparse kernel Jacobian against the dense one.

use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::gvalues::GValuesBoscolo;
use radkin::kernel::{SparseJacobian, Workspace, jacobian, jacobian_sparse};
use radkin::matrices::SimMatrices;
use radkin::physics::BeamCollection;
use serde_json::Value;

const CONFIGS: [&str; 3] = [
    "data/config_Labarbe.toml",
    "data/advanced_water_radiolysis.toml",
    "data/CONFIG_PAPER_01.toml",
];

fn matrices(path: &str) -> (SimMatrices, Vec<String>) {
    let cfg = Config::load(path).unwrap();
    let mut env = SimEnv::from_config(&cfg).unwrap();
    let beams = BeamCollection::from_config(&cfg.beam).unwrap();
    let (dose_rate, let_) = match beams.current() {
        Ok(b) => (Some(b.dose_rate_peak), beams.current_let().unwrap()),
        Err(_) => (Some(1e7), 0.201),
    };
    let gv = GValuesBoscolo::load(&radkin::consts::data_dir()).unwrap();
    env.apply_let(&gv, let_, dose_rate).unwrap();
    let labels = env.species.sorted_dyn_labels();
    (SimMatrices::build(&env).unwrap(), labels)
}

/// Densified `jacobian_sparse` and `jacobian` at `y`.
fn both(mat: &SimMatrices, sj: &SparseJacobian, y: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = y.len();
    let mut ws = Workspace::new(mat);
    let mut dense = vec![f64::NAN; n * n];
    jacobian(y, mat, &mut ws, &mut dense);
    let p = sj.pattern();
    let mut values = vec![f64::NAN; p.nnz()];
    jacobian_sparse(y, mat, sj, &mut ws, &mut values);
    let mut densified = vec![0.0; n * n];
    for j in 0..n {
        for k in p.col_ptr()[j]..p.col_ptr()[j + 1] {
            densified[p.row_idx()[k] + j * n] = values[k];
        }
    }
    (densified, dense)
}

/// The sparse fill does the same arithmetic as the dense one, in the same
/// order, so the entries agree exactly (the kernel golden tolerance is 1e-12).
fn assert_same(sparse: &[f64], dense: &[f64], what: &str) {
    for (i, (a, e)) in sparse.iter().zip(dense).enumerate() {
        assert!(a == e, "{what}[{i}]: {a} vs {e}");
    }
}

/// Deterministic pseudo-random states: concentrations over many decades,
/// a fraction of them negative (overshoot paths).
fn random_states(n: usize, count: usize, negative_share: f64) -> Vec<Vec<f64>> {
    let mut s = 0x9e3779b97f4a7c15u64;
    let mut next = move || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..count)
        .map(|_| {
            (0..n)
                .map(|_| {
                    let mag = 10f64.powf(-6.0 + 9.0 * next());
                    if next() < negative_share { -mag } else { mag }
                })
                .collect()
        })
        .collect()
}

#[test]
fn sparse_jacobian_matches_dense_on_random_states() {
    for path in CONFIGS {
        let (mat, _) = matrices(path);
        let sj = SparseJacobian::new(&mat);
        let n = mat.n_dyn;
        eprintln!(
            "{path}: n = {n}, nnz = {}, enzymatic = {}",
            sj.pattern().nnz(),
            mat.enzymatic.k_values.len()
        );
        assert!(
            sj.pattern().nnz() < n * n,
            "{path}: pattern should be sparse"
        );
        for (t, share) in [0.0, 0.3, 1.0].into_iter().enumerate() {
            for (k, y) in random_states(n, 20, share).into_iter().enumerate() {
                let (sparse, dense) = both(&mat, &sj, &y);
                assert_same(&sparse, &dense, &format!("{path} share {t} state {k}"));
            }
        }
    }
}

#[test]
fn sparse_jacobian_matches_dense_on_golden_states() {
    let g: Value =
        serde_json::from_str(&std::fs::read_to_string("tests/golden/derive.json").unwrap())
            .unwrap();
    for (name, path) in [("labarbe", CONFIGS[0]), ("pure_water", CONFIGS[1])] {
        let (mat, labels) = matrices(path);
        let sj = SparseJacobian::new(&mat);
        let vector = |state: &Value| {
            labels
                .iter()
                .map(|l| state[l].as_f64().unwrap())
                .collect::<Vec<_>>()
        };
        for y in [
            vector(&g[name]["states"]["synthetic"]),
            vector(&g[name]["overshoot_state"]),
        ] {
            let (sparse, dense) = both(&mat, &sj, &y);
            assert_same(&sparse, &dense, name);
        }
    }
}
