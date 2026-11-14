# Ticket 03: cli-resume-tests-docs

**Status:** ✅ Done

## Local Test Result

With CARGO_TARGET_DIR=.scratch/target: fmt --check clean, clippy -D warnings clean, `cargo test --release` all suites ok, no FAILED; tests/cli.rs `voxel_resume_from_checkpoint` ok (checkpoint written to --out, same-dir resume refused, resume past t_end refused, extension 1e-7 -> 2e-7 succeeds and yields end at 2e-7).

## Review Notes

Base 7f9ca51ccdce03784ba8a5048e075355fa30de30. Inline; no independent reviewer (CLI wiring, test and docs). `radkin voxel` now always writes `<out>/checkpoint.npz` and creates --out before the run; `--resume <ckpt>` reads it and refuses resuming into the same checkpoint path. The rtol equivalence, extension and grid-mismatch tests were landed in ticket 02 (tests/voxel_checkpoint.rs); this ticket adds the CLI test in tests/cli.rs instead of a new file. README documents checkpoint_every, --resume, the three uses, tolerance-level agreement, G-value baseline at the checkpoint time and np.load layout.

## Blockers / Challenges

None beyond the locked target/release/radkin.exe (another session's process), worked around with a scratch target dir.

## Commits

- 92415ae feat: radkin voxel --resume and checkpoint docs (ticket 03)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
