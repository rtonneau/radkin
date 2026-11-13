---
bump: minor
floor: minor
---
- `t_eval` points are now computed with scipy's LSODA dense output (order NQU, last-column rescale when the order drops) instead of DINTDY, matching radiopyo's interpolation formula. Values at `t_eval` times can shift by up to ~1e-7 relative.
- rsoda bumped to 312b742, which adds `Lsoda::dky` and `Lsoda::nordsieck`.
