# Ticket 03: kernel-optimization

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo clippy --release --all-targets --features timing -- -D warnings && cargo test --release` → exit 0, 57 tests passed (3 new kernel unit tests). Log `03-verify.log`.
`python -I tools/breakdown.py --repeat 10` → rhs share 21-35% (was 24-40%). Log `03-breakdown.json`.
A/B solve times, old vs new bench in alternation, 6x10 (`03-ab2.md`): every case faster, -4.1% to -9.7%, median about -8%. nfev and nst are identical on every case.

## Review Notes

Base: `6dc10f51eb73afb4f6b090e52fbea60b228a3bc0`.
- Kept: packed fixed-width mass-action reactions (`PackedK`: 3 slots padded with the factor-1 species, 4 rows padded with a dummy acc row of coef 0, ordered tails), and an RHS-only preparation (`fill_rhs`: `sw` only when there are negative products or a negative species, never `sl`). The products and per-species summation order are unchanged, so the results are bit-identical. Unit tests compare the bits against the reference composition (fill + compute_rates + species_rates + overshoot) on 200 states with up to 50% negative components, with one shared Workspace (stale-value check), and on a synthetic matrix with negative products. The tails are checked to be exercised.
- Dropped (recorded in docs/perf): flat CSR + fusion (slower), species-major gather (slower), unchecked indexing (only ~3%), reordered sums (cancellation vs the 1e-12 golden).
- The kernel stays allocation-free per call (packing happens in `Workspace::new`). The public API is unchanged (`reaction_rates`, `species_rates`, `Workspace::rates()` still filled).
- Jacobian (4-8%) left alone, as the ticket's ~15% threshold allows.

## Blockers / Challenges

The microbenchmarks are noisy (CPU clock): single runs varied 500-900 ns. I used old/new binaries run back to back with minimums. The first 3x7 A/B run had ±30% noise on 5 ms cases, and the 6x10 rerun was consistent.

## Commits

- 53a7ee9 perf(kernel): packed fixed-width reactions and RHS-only preparation in derive, bit-identical, ~8% faster solves (ticket 03)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
