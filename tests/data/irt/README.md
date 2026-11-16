# Reference data of the IRT stage

From Geant4-DNA 11.4.3 (`${DEV_DIR}/GEANT4/SIM/radkin`): 100 keV electrons through a
2 x 2 x 2 um water box in vacuum, physics option4, chemistry option3, IRT_syn.

- `g4_chemistry_dump.txt`: the reaction table Geant4 used at run time (species with
  their diffusion coefficient, charge and radius; reactions with their type, rate
  constants, radii and products). `tools/export_g4_chem.py` turns it into
  `data/chem_g4.toml`.
- `g4_events_molecules.csv`: the molecules of events 0 to 9 at 1.001 ps (positions in
  nm), from the pass that stops there.
- `g4_events_gt.csv`: the molecule counts of the same events in the G(t) pass at 1 ps
  and at every decade up to 1 us, under radkin's species labels.
- `g4_events_edep.csv`: the energy deposited by each event [eV].

The two passes share the physical stage of each event, not the random numbers of the
pre-chemical stage: their molecule counts at the start differ by a few molecules.
