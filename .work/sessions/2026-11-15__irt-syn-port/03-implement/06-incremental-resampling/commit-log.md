# Ticket 06: incremental-resampling

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_`: 43 passed, 0 failed (`irt_background` 5, `irt_incremental` 5, `irt_model` 7, `irt_sample` 9, `irt_schedule` 3, `irt_single_track` 8, `irt_walls` 6).

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure).

`cargo run --release --example irt_scaling -- --tracks 1,3,10,30,100` (and the larger cases run one by one):

```
tracks,molecules,resampling,wall_s,reactions,pair_times_drawn,e_aq_left
1,62,incremental,0.001,9,737,12
1,62,full,0.000,8,5048,13
3,252,incremental,0.001,44,5803,52
3,252,full,0.012,47,189700,47
10,869,incremental,0.006,162,55862,163
10,869,full,0.331,163,5655467,174
30,2607,incremental,0.045,489,396092,490
30,2607,full,7.747,509,119894211,490
100,8690,incremental,0.466,1839,4006507,1500
100,8690,full,280.558,1906,4133496418,1501
300,26070,incremental,4.451,6594,34093786,3872
1000,86900,incremental,61.940,31719,411565532,7945
```

Table and comments: `docs/perf/irt-scaling.md`. Raw files of the extra cases: `.scratch/tests/2026-11-15__irt-syn-port/06-scaling-*.csv`.

`python -I tools/compare_irt_geant4.py --events 100 --boundary outflow --irt-key 'resampling = "incremental"'`: RMS of z = 1.05 over 1065 values, every value of the decade table within 3 standard errors (H3O+ at 10 us: -2.3), 1.6 s for the 100 events (`.scratch/tests/2026-11-15__irt-syn-port/06-outflow-incremental/report.md`).

## Review Notes

Base: `05b4d9e19d90d71af2cd664fad32ac7e0d1f7e3e`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- `resampling = "incremental"` (`src/irt/incremental.rs`, `IncrementalStage`): reaction times are kept in a queue; a time is drawn for a pair when the later of its two particles appears (injected or created), and for a particle alone when it appears. After a reaction only the pairs of its products are drawn. `"full"` stays the default. `Engine` in `sim.rs` runs either rule behind one interface.
- Agreement with the full rule on the same seeds, every species at every record within 4 standard errors of the difference: three single tracks (largest |z| 3.16 over 270 values), five overlapping tracks with the pH 7 background (2.85), a 40 nm box with reflective walls (1.77).
- `examples/irt_scaling.rs` times both rules; the table is in `docs/perf/irt-scaling.md`.
- 100 overlapping tracks (8690 molecules) reach a 1 us handover in 0.47 s with the incremental rule, against 281 s with the full rule. 1000 tracks take 62 s.

How the incremental rule treats positions:
- A particle keeps the position of the last time it was looked at. It is moved to the present when a new particle looks for partners around it and when it reacts alone. The search allows for where it may have gone since (5 standard deviations of its displacement).
- A pair that fires is put at its reaction positions from the two positions its time was drawn from, over the time since then.
- Between reflective walls every pair is drawn again at the look-ahead limit of the box, as the full rule does at its box moves. Without it two particles that started further apart than the cutoff radius would never be paired.

Finding about Geant4 (noted in ticket 03): its IRT_syn keeps the reaction set between reaction steps and clears it only after a step limited by a process, so it is closer to this rule than to full resampling. The incremental rule matches the Geant4 reference at least as well as the full rule (RMS of z 1.05 against 1.00; OH at 10 us -1.2 against -2.3 standard errors).

Not in the ticket's file list: `src/irt/mod.rs`, `src/irt/sim.rs`, `src/irt/output.rs` (the engine switch and the rule name in `irt.json`), `docs/perf/README.md` (link).

## Blockers / Challenges

None.

## Commits

- 1d88b4c feat: IRT incremental resampling and scaling measure (ticket 06)

## Time Spent

13m (ticket-start.js to ticket-complete.js)
