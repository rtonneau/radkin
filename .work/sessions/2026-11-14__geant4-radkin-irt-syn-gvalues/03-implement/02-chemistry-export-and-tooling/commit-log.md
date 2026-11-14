# Ticket 02: chemistry-export-and-tooling

**Status:** ✅ Done

## Local Test Result

`python -I tools/run_sweep.py --only smoke`: 4 cases (event 0, 1 µm, {g4, air} × {10, 50 ns}), all exit 0 in 10–14 s.
- On the smoke seeds, G at the 10 ns handover equals the seed G (e⁻aq 2.989, •OH 2.988).
- `plot_gt.py` runs. IRT_syn at the handover vs the old plain IRT (pooled G(e⁻aq)): 2.974 vs 2.983 at 10 ns, 2.590 vs 2.602 at 50 ns.

## Review Notes

Base: 2112de0. Study files only (`.scratch`), no radkin code changes.
- `export_g4_chem.py` turns the dump into `configs/chem_g4.toml` plus `chem_g4.diffusion.json`. It gives 61 reactions, with reaction strings and rates identical (1e-5) to last session's hand-written `chem_option3.toml`.
- The 12 bulk H3Op(B)/OHm(B) reactions are dropped, each after checking that its track twin exists. H2O(B) reactions become k1/[H2O], NoneM ones first order. Autoionisation is added for pH 7.
- `make_configs.py`: seed time = handover, writes `[voxel.diffusion]` from the JSON. `run_sweep.py`: cases e<ev>_<chem>_<t>ns_dx<nm>. `plot_gt.py`: pools sum(G·E)/sum(E) from the net G columns.
- The aerated config is copied from the previous session.

## Blockers / Challenges

- An escaping slip in a generated Python line, fixed.

## Commits

- 2112de0 chore(gps): ticket 01-geant4-irt-syn-project done

## Time Spent

1m (ticket-start.js to ticket-complete.js)
