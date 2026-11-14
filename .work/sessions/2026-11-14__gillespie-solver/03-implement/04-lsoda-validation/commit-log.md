# Ticket 04: lsoda-validation

**Status:** ✅ Done

## Local Test Result

`cargo test --release --test ssa_validation -- --nocapture`:

```
reaction-only: 180 checks (78 on species above 100 molecules), worst |mean - lsoda| = 2.02 SEM, 79734 events per trajectory
continuous beam: 180 checks (76 on species above 100 molecules), worst |mean - lsoda| = 2.80 SEM, 105583 events per trajectory
pulsed beam: 180 checks (78 on species above 100 molecules), worst |mean - lsoda| = 2.59 SEM, 164860 events per trajectory
test result: ok. 3 passed; 0 failed; finished in 0.50s
```

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release` all passed (count in the completion output).

## Review Notes

Base: 9c211ebae2402f23b7042aeb0438f8f17942d632. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- Three cases on `data/advanced_water_radiolysis.toml`: reaction-only (radical initial state, `use_beam = false`), continuous beam (`run`), pulsed beam (`chunked_run`). The SSA is evaluated at 12 of the LSODA output times after 1 µs, so no interpolation is involved.
- Criterion `|mean - lsoda| <= 4 SEM_eff + 1 molecule`, `SEM_eff = max(sample SEM, sqrt(lsoda / N))`. No golden tolerance touched.
- Wrong-factor check by hand: with the propensity of `2 A` reactants multiplied by 2 in `SsaModel::propensity`, all three tests fail (H2O2: 4600 vs 4352 molecules at SEM 8.2 in the reaction-only case, i.e. 30 SEM). The change was reverted (`git diff src/ssa/model.rs` empty) and the tests pass again.
- Runtime: 0.5 s for the three cases.
- `docs/ssa.md` records the model, the cases, the criterion and the measured agreement.

Finding recorded in `docs/ssa.md`: with 2000-4000 trajectories the mean of H3O+ and OH- is about 4 molecules above LSODA in the beam cases, the same number of molecules at 1e-15 L and 1e-14 L. It is a real finite-volume effect (independent G-value sources make the net charge fluctuate; the water equilibrium is convex in it), not a propensity error: the reaction-only case shows none. The validation volume (1e-14 L) keeps it below half a SEM, so the 1-molecule bias allowance is not what absorbs it; the 4 SEM band is.

## Blockers / Challenges

None. The first exploration at 1e-15 L and 4000 trajectories showed deviations up to 11 SEM on H3O+/OH-; they were traced to the finite-volume effect above before the volume and criterion were fixed.

## Commits

- e3c5536 test: validate the SSA ensemble mean against LSODA at large volume (ticket 04)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
