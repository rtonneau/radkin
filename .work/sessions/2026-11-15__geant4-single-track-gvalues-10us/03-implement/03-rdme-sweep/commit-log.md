# Ticket 03: rdme-sweep

**Status:** ✅ Done

## Local Test Result

`python -I .scratch/tests/2026-11-15__geant4-single-track-gvalues-10us/tools/run_sweep.py --status`

```text
g4 1ps dx20 / g4 {1ps,100ps,1ns,10ns,100ns} dx10 / g4air ... dx10 o2L=2 / p01air ... dx10 o2L=2 / g4 1ps dx5: 10/10 each
ALL DONE
```

Wall times (runs/sweep.out): 20 nm 33 s per case. 10 nm: g4 about 235 s per case, aerated 545-580 s (2 jobs x 4 threads). 5 nm: 780-815 s per event (4 trajectories in sequence, 1 thread). Total about 4.2 h.

## Review Notes

Base bf006e8. Inline. Every output is in the gitignored scratch dir, so there is no repo commit.
- run_sweep.py runs every (case, event), skips finished runs (yields.csv plus mem.txt exit=0), and logs wall times. It was launched detached with Start-Process (pid 17116).
- 17 cases x 10 events x 4 trajectories, all exit=0, no reduction from the plan.
- Seedless background runs (`runs/_bg_rdme/`) were added for the analysis; see ticket 04.

## Blockers / Challenges

- None in the runs. The tool wait expired twice (2 h cap) and was re-armed. The analysis found a background drift (ticket 04).

## Commits

- 004414b chore(gps): ticket 02-configs-and-pilot done

## Time Spent

4h 09m (ticket-start.js to ticket-complete.js)
