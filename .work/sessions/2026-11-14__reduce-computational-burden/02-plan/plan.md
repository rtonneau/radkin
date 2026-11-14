# Implementation Plan

**Session:** reduce-computational-burden
**Date:** 2026-11-14T13:15:30.411Z
**Estimated effort:** 1-2 days

## Strategy

Measure first, then apply the three approved levers one at a time, each measured against the previous kept state with the 5% keep rule and the accuracy gate. Benchmark and accuracy harness come first so every lever has numbers.

## Tickets Overview

- **Ticket 1:** Voxel benchmark (`examples/voxel_bench.rs`) on the Geant4-seeded case plus the tight-rtol accuracy gate and the baseline in `docs/perf/voxel.md`.
- **Ticket 2:** Lean per-voxel solve: per-thread reuse of solver state, no t_eval output, carried step size.
- **Ticket 3:** Fuse adjacent half-reactions in the dt/2 path of step-doubling.
- **Ticket 4:** Embedded splitting error estimate replacing step-doubling, with an ADR.

## Sequencing Rationale

The benchmark gates every lever. The lean solve is independent of the driver and makes later measurements less noisy. Fusion is a local driver change; the embedded estimator rewrites the same driver loop and may subsume fusion, so it comes last.

## Risks & Mitigation

- **Risk:** carrying step size or reusing solver state changes results or breaks determinism across thread counts → **Mitigation:** state is per voxel (not per thread) for the carried step; parallel vs serial bit-identity test kept; accuracy gate.
- **Risk:** embedded estimator under-estimates error and accepts too-large steps → **Mitigation:** accuracy gate vs tight reference; existing `error_decreases_with_tolerance` test; dropped if gate fails.
- **Risk:** lean solve needs rsoda API → **Mitigation:** change rsoda in its repo, push, bump `rev`; temporary `[patch]` only locally.
- **Risk:** checkpoints must stay compatible → **Mitigation:** carried step size is not persisted (resume re-estimates), format unchanged; `tests/voxel_checkpoint.rs` green.

## Assumptions

- `data/voxel_example.toml` with `data/voxel_example_molecules.csv` (Geant4 molecules) is the Geant4-seeded benchmark case; a longer end time may be set in the bench to make it representative.
- Phase-1 (unit-cell) paths and radiopyo golden tests are not affected.
