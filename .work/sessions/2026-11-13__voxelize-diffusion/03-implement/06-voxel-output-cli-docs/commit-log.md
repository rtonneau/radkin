# Ticket 06: voxel-output-cli-docs

**Status:** ✅ Done

## Local Test Result

fmt --check clean; clippy -D warnings clean with and without `parallel`; `cargo test --release` 15 suites ok (cli 4 passed incl. voxel example). `radkin voxel data/voxel_example.toml` → 768 voxels, 57 steps, 0 rejected, 1.0 s, 4 snapshots. `python -I tools/compare_geant4.py` renders G(t) with a dummy Geant4 CSV.

## Review Notes

Base bea65f6. `VoxelSimulation` (chemistry via Config::load + `[voxel]` + `[voxel.run]`), `VoxelRun::track_totals` records domain sums after each accepted step and injection. Outputs: yields.csv (N_<sp>, G_<sp> = ΔN/(E_injected/100 eV)), .npy snapshots in mol/L with C-order shape (nz,ny,nx) (header 64-byte aligned, tested), snapshots.json. CLI `radkin voxel CONFIG --out DIR`. Example config plus a synthetic 66-molecule track (marked illustrative). README Phase-2 section with the mean-field / voxel-size caveat; AGENTS.md layout updated. Self-review of the plot found the pre-injection row at the injection time getting G=0; fixed by keeping only the post-injection row for duplicate times (tested).

## Blockers / Challenges

Pre-existing, not touched: output::RSODA_REV says 312b742 while Cargo.toml pins 2f47645.

## Commits

- ede169b feat(voxel): radkin voxel CLI, yields/G(t) CSV, npy snapshots, example, Geant4 comparison tool and docs (ticket 06)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
