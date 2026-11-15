# Ticket 02: configs-and-pilot

**Status:** ✅ Done

## Local Test Result

`cargo build --release && python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/make_configs.py --check`

```text
Finished `release` profile [optimized]
g4: ok
g4air: ok
p01air: ok
OK
```

Pilot results are in `.scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/pilot.md`. Event 0, 1 ps, 10 nm: 36 s per trajectory, and G at 1.001 ps equals the seed G (e_aq 3.86). The yields cover 1 ps to 10 us.

## Review Notes

Base bf006e8 (tree clean apart from the session record). Inline, no reviewer: every file is in the gitignored scratch dir, so nothing goes into the repo.
- `tools/make_configs.py`: per (chem, handover, dx, event) config with `solver = "ssa"`, 4 trajectories, seed 1000 + event, grid [-1000, 1000] nm^3, t_start = handover, t_end = 10 us, energy_ev = edep. Seeds are clamped into the box (`seeds/<pass>.csv`): IRT_syn keeps molecules on the closed box, while the grid's upper faces are open.
- Chemistries: chem_g4 (copied from 11-14), chem_g4_air (the same plus O2 = 52e-6 M), chem_paper01_aerated.
- O2 level 2 for both aerated chemistries: it matches level 0 within 1.5 SEM, with 13.7x fewer hops (pilot.md). radkin's own suggestion is 0.
- Sweep sizing: 10/20 nm at 2 jobs x 4 threads; 5 nm at 1 job x 1 thread (21.5 GB per trajectory).

## Blockers / Challenges

- The first 5 nm pilot failed with an OOM (15.36 GB allocation) while it ran next to the other pilots. Alone it passes at about 21.5 GB. Part of that memory is the dense f64 seed field that the voxel loader keeps alongside the RDME counts. A sparse seed representation would halve the 5 nm memory; this is noted as an idea for later and is out of scope here.
- The working set under-reports committed memory on Windows. The PeakPagedMemorySize64 commit probe in runmem.ps1 reported 0, so the 21.5 GB figure was read from PrivateMemorySize64 during the run.

## Commits

- b46df9f chore(gps): ticket 01-geant4-segment-passes done

## Time Spent

35m (ticket-start.js to ticket-complete.js)
