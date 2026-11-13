# Ticket 06: pin-rsoda-and-final-report

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo clippy --release --all-targets --features timing -- -D warnings && cargo test --release` → exit 0, 57 passed (`06-verify.log`). `git grep -n "patch\|path = " Cargo.toml` → only `path = "src/bin/radkin.rs"` (the bin target), no rsoda patch or path. Cargo.lock: rsoda `git+...?rev=f5d1312`.
Final A/B, baseline vs final, 6x10 (`06-ab-final.md`): -16.5% to -25.4% on every case, nfev and nst identical.

## Review Notes

Base: `9fd3d1e6c1e3ec34c87134a78f289535788d2341`.
- rsoda: version 0.2.1 (`f5d1312`) on `perf/per-step-overhead`, pushed. PR https://github.com/rtonneau/rsoda/pull/4 is open and not merged (merging is the user's call). The README already documents the timing feature (ticket 1), and rsoda has no changelog file.
- radkin: Cargo.toml pins rev `f5d1312`, and the temporary `[patch]` is removed. docs/perf has the Final section: the A/B vs baseline, the final breakdown, and what is left and why. It links the idea through the session note at the top.

## Blockers / Challenges

None. rsoda PR CI had not been checked when this was written. It must be green before merge.

## Commits

- d04e725 chore(deps): pin rsoda f5d1312 (0.2.1); final perf report, 16-25% faster solves (ticket 06)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
