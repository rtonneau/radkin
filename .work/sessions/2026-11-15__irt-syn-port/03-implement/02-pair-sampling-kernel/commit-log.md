# Ticket 02: pair-sampling-kernel

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_sample`: 9 passed, 0 failed (`tests/irt_sample.rs`). `cargo test --release --lib rng`: 10 passed (3 new: `det_exp`, `norm_quantile`, `normal`).

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure).

## Review Notes

Base: `7ca6ba299997557d203f4d79c652b1a031114f20`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- `pair_time` follows `G4DiffusionControlledReactionModel::GetTimeToEncounter`: first-encounter time from `W = Reff / r`, then for a partially diffusion-controlled pair the acceptance `k_act / (k_act + k_dif)` and the second time.
- `contact_reacts` gives the contact probability of a pair found inside `Reff` at the start.
- `reaction_positions` and `product_position` follow `G4DNAMakeReaction`, with two deliberate differences where Geant4 departs from the formulas it implements: the pair coordinate diffuses with the standard deviation `sqrt(sigma)` (Geant4 passes the variance, which in its units is about zero), and the angle of the separation at the reaction is taken from the initial separation (Geant4 takes it from the z axis). Both are stated in the module doc.
- `erfc_inv` is Geant4's (`-NormQuantile(x / 2) / sqrt(2)`, Wichura AS 241); `erfc(erfc_inv(x))` returns `x` within 1e-12 relative from 1e-300 to 1, and 1e-12 absolute up to 2.
- Empirical reaction probabilities (400 000 pairs, 4 standard errors) match the analytic ones at 5 to 6 times from 1 ps to 1 us and in the long-time limit: fully diffusion-controlled (Smoluchowski), partially diffusion-controlled and ion pairs (radiation boundary). The second time of Geant4 is exactly the radiation-boundary law, not only its long-time limit.
- The kernel takes `&mut SsaRng` and plain values, allocates nothing, and uses IEEE basic operations, `sqrt`, `det_ln` and the new `det_exp` only.

Added to `src/ssa/rng.rs`: `det_exp`, `norm_quantile`, `SsaRng::uniform_inside`, `SsaRng::normal` (one uniform per normal, by inversion). `effective_radius` of ticket 1 now uses `det_exp`, so the radii do not depend on the platform libm; the 52 radii still match the Geant4 dump.

Follow-up to ticket 1 in this commit: the test functions of `tests/irt_model.rs` are renamed `irt_model_*` (and those of this ticket `irt_sample_*`). `cargo test --release irt_model` filters on test names, so the Verification Step of ticket 1 as written ran no test; I had run it with `--test irt_model`. With the prefix both commands run the 7 tests.

## Blockers / Challenges

- The analytic radiation-boundary reference of the test overflowed at long times (`exp(...) * erfc(...)`); it now uses the scaled form `exp(-x^2) erfcx(x + y)`.
- A reference quantile typed from memory was wrong (the code was right, checked with mpmath).
- Clippy rejects the published 20-digit AS 241 coefficients (`excessive_precision`); the lint is allowed on `norm_quantile` and `det_exp` so the published digits stay.

## Commits

- f5d86ef feat: IRT pair sampling kernel (ticket 02)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
