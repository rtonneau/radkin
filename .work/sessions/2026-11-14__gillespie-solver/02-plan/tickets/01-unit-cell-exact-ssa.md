# Ticket 01: unit-cell-exact-ssa

**Acceptance Criteria:**
- [ ] `src/ssa/` exists with a propensity model built from `SimMatrices` and a volume [L]: mass-action reactions use `k · Π[fixed] · N_A·V · Π falling-factorial(dynamic counts) / (N_A·V)^m`, Michaelis–Menten reactions use their rate law, fixed species keep constant concentrations
- [ ] A reaction can fire only when every species it consumes on balance (negative products included) is present in sufficient number
- [ ] A + A (written `2 A -> ...` or `A + A -> ...`) has propensity `k/(N_A·V)·n·(n−1)` and removes 2 molecules per firing (unit test)
- [ ] Deterministic RNG: ChaCha8, stream index = trajectory index; uniform and exponential samplers use no libm call
- [ ] One exact-SSA trajectory (direct method) runs on a reaction-only case from given initial counts, records counts at requested output times, and is identical for the same seed and stream
- [ ] `cargo fmt --check`, clippy `-D warnings` and `cargo test --release` pass

**Files to Touch:**
- `Cargo.toml`
- `Cargo.lock`
- `src/lib.rs`
- `src/ssa/mod.rs`
- `src/ssa/rng.rs`
- `src/ssa/model.rs`
- `src/ssa/exact.rs`
- `tests/ssa_exact.rs`

**Verification Step:**

Run:
```bash
cargo test --release ssa
```

Expected:
All `ssa` unit and integration tests pass (0 failed).

**Notes:**

`mat.k.reactant_index` lists one slot per unit of stoichiometry (dynamic species `< n_dyn`, fixed species above); `mat.k.stoichiometry` holds the net integer change. State counts are per dynamic species in kernel order. Initial counts for the unit cell: floor of `c·N_A·V` plus a Bernoulli draw of the fraction.
