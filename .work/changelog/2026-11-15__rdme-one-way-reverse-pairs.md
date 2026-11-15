---
bump: minor
floor: minor
---
- The RDME mesh correction also keeps the equilibrium of reversible pairs written as two one-way reactions (e.g. `H3O_plus + OH_minus -> H2O` and `H2O -> H3O_plus + OH_minus`, as in Geant4 exports). Exact reverses are matched on their dynamic species and logged as inferred. Ambiguous candidates get a warning and are left unpaired.
