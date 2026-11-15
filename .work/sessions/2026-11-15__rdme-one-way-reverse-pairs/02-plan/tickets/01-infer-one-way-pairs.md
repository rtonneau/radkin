# Ticket 01: infer-one-way-pairs

**Model:** heavy
**Effort:** medium

**Acceptance Criteria:**
- [ ] `SimMatrices` adds inferred pairs after the declared ones: two mass-action one-way reactions A, B with dynamic-species multisets (fixed species dropped) reactants(A) = products(B) and products(A) = reactants(B). The pair is formed only when each has exactly one such candidate, and each pair is marked as inferred.
- [ ] A reaction with more than one candidate is not paired, and a warning naming the candidates is logged at RDME setup.
- [ ] One-way reactions are never paired with a direction of a declared `<->` pair: in CONFIG_PAPER_01 the existing test still sees exactly the 9 declared pairs, and the one-way `e_aq + H2O -> H_r + OH_minus` stays unpaired.
- [ ] A unit test on a `chem_g4.toml`-style network (`H3O_plus + OH_minus -> H2O`, `H2O -> H3O_plus + OH_minus`) infers the pair, and `RdmeSetup` gives both directions the shared factor.
- [ ] `log_reversible_pairs` labels inferred pairs "inferred from one-way reactions".
- [ ] ADR 0006 body, `docs/ssa.md` and the glossary "Reversible pair" entry describe the inference. fmt, clippy -D warnings and `cargo test --release` pass.

**Files to Touch:**
- `src/matrices.rs`
- `src/ssa/rdme.rs`
- `src/ssa/voxel_sim.rs`
- `docs/ssa.md`
- `.work/adr/0006-rdme-reversible-pair-correction.md`
- `.work/GLOSSARY.md`

**Verification Step:**

Run:
```bash
cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
```

Expected:
All tests pass; clippy clean.

**Notes:**

The user chose matching on dynamic species only, with no config key (always on under hellander), and unit tests only (no reruns). See `.work/sessions/2026-11-15__rdme-one-way-reverse-pairs/01-grill/resume.md`. The study chemistry is `.scratch/tests/2026-11-15__rdme-reversible-mesh-correction/configs/chem_g4.toml`; copy the minimal network into a test fixture rather than depending on .scratch.
