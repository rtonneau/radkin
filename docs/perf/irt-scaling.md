# IRT stage: run time against the number of overlapping tracks

Measured on 2026-11-15 (gps session `2026-11-15__irt-syn-port`, ticket 06) with
`examples/irt_scaling.rs`:

```text
cargo run --release --example irt_scaling -- --tracks 1,3,10,30,100
```

## Case

- Tracks: the ten reference events of `tests/data/irt/` (2 um segments of 100 keV
  electron tracks, 87 molecules each on average at 1.001 ps), used in turn.
- All tracks are injected together and moved sideways by a random offset inside a
  square of 1000 nm: they overlap.
- Chemistry `data/chem_g4.toml` with its background reactions (water, pH 7), free
  space, `start_contact = false`.
- One trajectory per case, one thread, wall time of the run without the loading.
  16-core desktop, Windows 11.

## To a 1 us handover

| tracks | molecules | incremental | full | reactions |
|---|---|---|---|---|
| 1 | 62 | 1 ms | below 1 ms | 9 |
| 3 | 252 | 1 ms | 12 ms | 44 |
| 10 | 869 | 6 ms | 0.33 s | 162 |
| 30 | 2 607 | 47 ms | 7.7 s | 489 |
| 100 | 8 690 | 0.47 s | 281 s | 1 839 |
| 300 | 26 070 | 4.5 s | not run | 6 594 |
| 1000 | 86 900 | 62 s | not run | 31 719 |

- The target of the session, about 100 overlapping tracks to the handover in
  minutes, is met by the incremental rule in half a second per trajectory.
- Full resampling draws every pair again after each reaction: its time grows about
  as the cube of the number of molecules (4.1e9 reaction times drawn at 100 tracks,
  against 4.0e6). It is a cross-check (`docs/irt.md`) and stays usable up to about
  30 tracks.
- Incremental resampling grows about as the square: every new molecule looks at
  every other one, and the first injection pairs every molecule with every other.
  At 1000 tracks the run spends most of its time there.
- The two rules give close chemistry: 1500 and 1501 e_aq left at 100 tracks (one
  trajectory each). The incremental rule is the default; the full rule is a few per
  cent off Geant4 on H2O2 (`docs/irt.md`).

## Other cases, incremental rule

| case | tracks | molecules | time | reactions |
|---|---|---|---|---|
| to 10 us | 10 | 869 | 16 ms | 195 |
| to 10 us | 100 | 8 690 | 1.8 s | 3 052 |
| to 10 us | 300 | 26 070 | 19 s | 12 125 |
| to 1 us, square of 200 nm | 100 | 8 690 | 1.1 s | 3 062 |

- A later handover costs more than its extra reactions: the cutoff radius of the
  pair search grows with the time left (1.4 um for 10 us, against 0.44 um for 1 us),
  so each new molecule draws more reaction times.
- A denser bundle reacts more (3062 reactions against 1839) and costs twice as much.

## Where the time goes, and what would come next

- Incremental: the search of the partners of a new molecule scans every molecule
  (a distance test each), then moves the ones in reach to the present and draws
  their reaction times. A persistent grid of cells would make the scan local; it
  is not needed below about 300 tracks.
- Full: the reaction times of all pairs inside the cutoff radius, at every
  reaction. Nothing short of keeping the times (the incremental rule) helps.
- Trajectories run in parallel (one per thread); the numbers above are for one.
