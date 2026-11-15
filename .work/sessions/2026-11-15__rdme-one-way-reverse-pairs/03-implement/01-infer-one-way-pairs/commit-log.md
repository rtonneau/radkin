# Ticket 01: infer-one-way-pairs

**Status:** ✅ Done

## Local Test Result

`cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release`: exit 0. fmt clean, clippy clean, every test binary `ok` (lib: 151 passed, 0 failed). Full log: `.scratch/tests/2026-11-15__rdme-one-way-reverse-pairs/evidence/01-verification.log`. A later edit touched only `docs/ssa.md` prose.

## Review Notes

- AC1: `SimMatrices::reversible_pairs` is now `Vec<ReactionPair { forward, backward, inferred }>`; `infer_one_way_pairs` appends inferred pairs after the declared ones. Candidates are one-way mass-action reactions (not a direction of a declared pair, no negative product, non-zero net change) with sorted dynamic reactant/product multisets; a pair forms only when each is the single candidate of the other. Tests `one_way_reverses_are_inferred_as_pairs`, `inferred_pairs_match_dynamic_species_only`.
- AC2: `SimMatrices::ambiguous_reverses` lists each reaction with more than one candidate and its candidates; `log_reversible_pairs` (voxel_sim, under hellander) warns naming them. Test `a_reaction_with_several_reverses_is_not_paired`.
- AC3: one-way reactions are never matched with declared directions (test `one_way_reactions_are_not_paired_with_a_declared_direction`); in CONFIG_PAPER_01 the 9 declared pairs are still the first 9, not inferred, and `e_aq + H2O -> H_r + OH_minus` stays unpaired. Deviation: CONFIG_PAPER_01 also holds two genuine one-way reverse pairs (`H_r + H2O -> H2 + OH_r` / `OH_r + H2 -> H_r + H2O`, `O_minus + O2 -> O3_minus` / `O3_minus -> O_minus + O2`), which AC1's rule infers, so the total is 11, not 9. The CONFIG_PAPER_01 tests in `matrices.rs` and `ssa/model.rs` now assert the 9 declared plus these 2 inferred.
- AC4: `rdme.rs` test `inferred_pairs_share_the_factor_of_their_bimolecular_direction` on a chem_g4-style network (inline fixture, no .scratch dependency) checks `ReversiblePair::Shared { inferred: true }` and the shared factor on levels 0 and 1.
- AC5: `ReversiblePair` variants carry `inferred`; `log_reversible_pairs` names an inferred pair by both reactions with "(inferred from one-way reactions)".
- AC6: ADR 0006 body, `docs/ssa.md` (Reversible pairs section) and the glossary entry updated.

## Blockers / Challenges

The AC3 wording "exactly the 9 declared pairs" conflicts with AC1 on CONFIG_PAPER_01, which has two inferable one-way pairs (see Review Notes). AC1's rule was kept and the tests assert 9 declared + 2 inferred; the parent should confirm this reading. Voxel hellander runs of CONFIG_PAPER_01 therefore now also share factors for those two pairs.

### Parent review and acceptance

- Base 9aa20558cb9deb47ffc4588e5c737ac6f9352554; pre-existing changes: gps session state only. Correction passes used: 1 (parent).
- AC3 reading (parent): its intent is that one-way reactions are never paired with a declared direction. "Exactly 9" was a wrong expectation: CONFIG_PAPER_01 has 2 real one-way reverse pairs (`H_r + H2O` / `OH_r + H2`, `O_minus + O2` / `O3_minus`), which the user's dynamic-species rule pairs. Their voxel hellander results change accordingly.
- Independent review (base 9aa2055): P2 `src/ssa/voxel_sim.rs` three log literals lost their `\` continuations; fixed by the parent. Re-review: no findings.
- Verification after the fix: fmt, clippy and cargo test --release pass (`evidence/01-parent-fix.log`).
- Accepted with known gaps: no tests on the log text, and no RdmeSetup test for an ambiguous case (the matrices-level ambiguity test exists). Requested model opus/medium.

## Commits

- dae9ecc feat: infer reversible pairs from exact one-way reverses (ticket 01)

## Time Spent

12m (ticket-start.js to ticket-complete.js)
