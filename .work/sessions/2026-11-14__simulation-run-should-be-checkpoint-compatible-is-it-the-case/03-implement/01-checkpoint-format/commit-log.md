# Ticket 01: checkpoint-format

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check` clean; `cargo clippy --release --all-targets -- -D warnings` clean; `cargo test --release --lib`: 52 passed (4 checkpoint tests). Manual: `np.load` on a written checkpoint lists state (2,2,2,3), time, grid, seeds, background, meta.json.

## Review Notes

Base f944236a30da582f38988d1da87a35fb01b75aa6; pre-existing changes: none in src (other session's .work files untouched). Inline, no independent reviewer (localized new module, checked directly). All 4 criteria met. Choices: floats stored in .npy entries (t, dt, grid, seed times) for bit-exact round trip without changing serde_json float parsing; `check_compatible` also compares the grid origin (it moves seed placement); new `RadkinError::Checkpoint`; dependency `zip` 2 without default features (stored entries only, pure Rust).

## Blockers / Challenges

A running radkin.exe (PID 9292, another session's sweep) locks target/release/radkin.exe, so `cargo test` fails to relink the bin; lib tests run with `--lib`. Not killed.

## Commits

- a645a3f feat: Voxel checkpoint .npz format (ticket 01)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
