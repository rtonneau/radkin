# RDME mesh correction of reversible pairs

_ADR 0006, recorded 2026-11-15 (session 2026-11-15__rdme-reversible-mesh-correction)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

Hellander mesh correction (ADR 0004) scaled only bimolecular propensities, which broke detailed balance for `<->` pairs such as water autoionisation: the pH 7 background drifted with voxel size. Both directions of each declared reversible pair now share the factor of their bimolecular direction, so K_eq is kept (Fange 2010; Hellander 2012). Pairs that are bimolecular on both sides keep separate factors and log the K_eq shift, and nothing is inferred from one-way reactions.
