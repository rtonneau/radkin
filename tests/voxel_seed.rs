use radkin::chemistry::env::SimEnv;
use radkin::config::Config;
use radkin::consts::data_dir;
use radkin::voxel::{VoxelConfig, load_seeds};

#[test]
fn seed_files_load_through_the_config() {
    let dir = std::env::temp_dir().join(format!("radkin-voxel-seed-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("track.csv"),
        "species,x,y,z,event_id\ne_aq^-1,1,1,1,0\n°OH^0,2,2,2,0\nH3O^1,3,3,3,1\nH2O^0,0,0,0,1\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("run.toml"),
        r#"
        [voxel.grid]
        dims = [4, 4, 4]
        spacing_nm = 1.0
        [voxel.geant4]
        ignore = ["H2O^0"]
        [[voxel.seed]]
        file = "track.csv"
        time = 1e-12
        events = [0]
        [[voxel.seed]]
        file = "track.csv"
        time = 1e-9
        "#,
    )
    .unwrap();
    let vc = VoxelConfig::from_toml_file(dir.join("run.toml")).unwrap();
    let cfg = Config::load(data_dir().join("CONFIG_PAPER_01.toml")).unwrap();
    let env = SimEnv::from_config(&cfg).unwrap();
    let seeds = load_seeds(&vc.seeds, &vc.grid, &env.species, &vc.policy).unwrap();
    assert_eq!(seeds.len(), 2);
    assert_eq!(seeds[0].report.placed, 2);
    assert_eq!(seeds[1].report.placed, 3);
    assert_eq!(seeds[1].report.ignored, 1);
    let c = vc.grid.molecule_concentration() * 1e6;
    let total: f64 = seeds[1].amount.data.iter().sum();
    assert!((total / c - 3.0).abs() < 1e-12);
    std::fs::remove_dir_all(&dir).ok();
}
