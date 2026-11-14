---
bump: patch
floor: patch
---
- `radkin voxel` runs faster: about 21% at 16 threads and 26% single-threaded on the Geant4-seeded example, with unchanged accuracy (leaner per-voxel reaction solves and fused half-reactions in the step-doubling control).
- New `examples/voxel_bench.rs` voxel benchmark with an accuracy check against a tight-tolerance reference; results in `docs/perf/voxel.md`.
