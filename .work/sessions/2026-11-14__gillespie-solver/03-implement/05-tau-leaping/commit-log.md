# Ticket 05: tau-leaping

**Status:** ✅ Done

## Local Test Result

`cargo test --release ssa`: all `ssa` unit tests (rng with the Poisson sampler, model, schedule, config, ensemble, tau) and integration tests (`ssa_exact`, `ssa_beam`, `ssa_cli`, `ssa_validation`, `ssa_tau`) pass.

`cargo test --release --test ssa_tau -- --nocapture`:

```
epsilon 0.001: worst mean difference 2.38 combined SEM, worst std difference 40.9 %, events ratio 1.0000
epsilon 0.03: worst mean difference 16.22 combined SEM, worst std difference 29.1 %, events ratio 1.0002
pulsed water: worst mean difference 3.75 combined SEM, worst std difference 43.4 %, events ratio 1.0047
test result: ok. 3 passed; 0 failed
```

Full checks: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release`: 177 passed, 0 failed.

## Review Notes

Base: 2c447eb3b183486c0c5b970d4ceb3b9a385c4638. Pre-existing changes: none outside `.work/`.

Checked against the Acceptance Criteria:
- `src/ssa/tau.rs`: CGP 2006 step selection (epsilon, g_i from the highest-order role of each species), critical reactions (fewer than 10 firings left) fired one at a time, fallback to 100 exact steps when the leap is below 10/a0.
- Counts never go negative (a trial state is drawn first; a negative count rejects the leap and halves it); the step is capped at the end of the span, so beam transitions and output times are not crossed. Unit test over 50 streams checks non-negativity, three conservation laws and the exact end time; the integration test checks that nothing is produced between pulses.
- Poisson sampler without libm (`SsaRng::poisson`: exponential arrivals below 10, Hormann PTRS above, `ln_factorial` from `det_ln`), tested on mean, variance and a chi-square against the exact probabilities on both sides of the switch.
- Speed-up measured with `examples/ssa_bench.rs` and recorded in `docs/ssa.md`, including the water case where the fast equilibria leave no gain.

Deviation from the ticket's fourth criterion (my decision, flagged for the user): the criterion said the tau-leaping and exact means agree within 4 combined SEM on a large-count case. That is not true of tau-leaping at its default accuracy: the mean has a systematic error proportional to epsilon (32 combined SEM at 1e-12 L with epsilon 0.03, falling to 11, 3.7 and 2.4 at 0.01, 0.003 and 0.001; 0.36 % of the species peak at 0.03). The test therefore checks (a) agreement within 4 combined SEM with epsilon = 0.001 and (b) agreement within 4 combined SEM plus epsilon times the species peak at the default 0.03. To make (a) possible and to give users the accuracy knob, `[ssa] tau_epsilon` (default 0.03) was added; it is not in the grill.

Also corrected in `docs/ssa.md`: the single-core rate of the exact SSA is about 5 M events/s on pure water (the earlier figure came from a 16-thread run).

## Blockers / Challenges

The first comparison showed tau-leaping 32 SEM away from the exact SSA at very large counts. It was checked to be the method's first-order error (it shrinks in proportion to epsilon) and not a bug, then the test criterion was written accordingly.

## Commits

- e96dbd5 feat: opt-in tau-leaping (Cao-Gillespie-Petzold) with a tau_epsilon key (ticket 05)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
