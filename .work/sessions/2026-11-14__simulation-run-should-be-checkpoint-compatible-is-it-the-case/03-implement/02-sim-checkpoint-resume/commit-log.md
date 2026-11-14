# Ticket 02: sim-checkpoint-resume

**Status:** ✅ Done

## Local Test Result

With CARGO_TARGET_DIR=.scratch/target: `cargo fmt --check` clean, `cargo clippy --release --all-targets -- -D warnings` clean, `cargo test --release` all suites ok (lib 53 incl. checkpoint_every_setting; tests/voxel_checkpoint.rs 3 passed: resume vs uninterrupted within rtol 1e-3/atol 1e-4 scaled diff <= 1, extend 1e-10 -> 1e-9, seed skipping, new seed after checkpoint, missing seed error, grid mismatch and t_end <= t_ckpt errors).

## Review Notes

Base b58329d7830de356143685584848d105e1e2b444. Inline; no independent reviewer (single module plus tests, checked directly against criteria). Criteria met: `checkpoint_every` (integer >= 1) plus final checkpoint whenever `RunControl::checkpoint` is set; `run_with(&RunControl{checkpoint, resume})` resumes t, state, dt (`VoxelRun::set_dt`), background, Rannacher pending; seeds matched by (time, fingerprint), unmatched earlier seed is a Checkpoint error, unmatched records warned; settings, chemistry hash (Debug of SimEnv, deterministic) and seed list compared with log::warn. Deviation: library integration tests added in tests/voxel_checkpoint.rs in this ticket instead of ticket 3 (ticket 3 keeps CLI test and docs).

## Blockers / Challenges

radkin.exe of another session still locks target/release; full test suite run in .scratch/target.

## Commits

- 9979269 feat: Checkpoint writing and resume of voxel runs (ticket 02)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
