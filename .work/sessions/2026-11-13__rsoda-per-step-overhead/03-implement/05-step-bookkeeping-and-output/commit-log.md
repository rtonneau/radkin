# Ticket 05: step-bookkeeping-and-output

**Status:** ✅ Done

## Local Test Result

rsoda (`aa5c294`, pushed to `perf/per-step-overhead`): `cargo fmt --check && cargo clippy --release --all-targets --all-features -- -D warnings && cargo test --release && cargo test --release --all-features && cargo doc --no-deps --all-features` → exit 0, 267 tests over both runs, golden included, plus the new `predict_and_retract_match_the_indexed_sweeps_bit_for_bit` (log `05-rsoda-verify.log`).
radkin with the local `[patch]` to rsoda (not committed): `cargo clippy --release --all-targets -- -D warnings && cargo test --release` → exit 0, 57 passed (`05-radkin-verify.log`).
A/B vs the ticket-3 bench, 6x10 (`05-ab-slices.md`): -5.8% to -18.4% on every case. nfev and nst identical everywhere.

## Review Notes

Base: radkin `06c54d20f30030039b5212ecf6bf6215c832e01d`, rsoda `08a2abf`.
- Profiling: a temporary rdtsc probe per DSTODA label plus the driver (patches kept in scratch `05-probe-*.diff`, not committed) found rsoda's own work at ~32% of its run time: predict 9%, corrector vector loops ~9%, YH update 7%, driver 3%. radkin's output storage (solve minus step time) is 5-9%.
- Kept (rsoda): pascal_sweeps on column slices for predict/retract, and zipped loops for the functional and chord corrector, the YH update and the L178 rescale. Same ops and order, so bit-identical: rsoda golden passes with exact counters, and the unit test compares the bits against the indexed sweeps for nyh 1/3/22 and nq 1..12. Step heuristics untouched (nst/nfe/nje identical).
- Dropped (below 5%): 4-lane vmnorm (-0.6..-3.6%) and flat output storage (0..-7.6%, mostly -1..-3.6%, and a public API change). Both recorded in docs/perf with numbers.
- radkin has no source change in this ticket (docs only). The rsoda rev pin is ticket 6.

## Blockers / Challenges

There is no sampling profiler (xperf needs elevation), so I used a scratch rdtsc probe. A first 'no storage' upper-bound hack broke chunk continuation (nfev changed on chunked cases), so only the run-mode rows of that experiment count, as noted in docs/perf.

## Commits

- c51ac6a docs(perf): rsoda step-loop slices kept (6-18% faster), vmnorm lanes and flat output dropped (ticket 05)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
