# Ticket 01: rate-correction-config

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] `SsaConfig` has `rate_correction: Option<RateCorrection>` (`RateCorrection::{Hellander, None}`); unset is `None` (meaning the default, Hellander, with fallbacks).
- [ ] `rate_correction = "hellander"` and `"none"` parse; any other value or type is a config error naming `ssa.rate_correction`; the unknown-key message lists it.
- [ ] A unit-cell run with `solver = "ssa"` and an explicit `rate_correction = "hellander"` fails with a config error; unset or `"none"` is accepted.
- [ ] The module doc of `src/ssa/config.rs` shows the key.

**Files to Touch:**
- `src/ssa/config.rs`
- `src/ssa/mod.rs`
- `src/bin/radkin.rs` (or the unit-cell ensemble entry, wherever `require_volume` is checked)

**Verification Step:**

Run:
```bash
cargo test --release ssa::config
```

Expected:
All config tests pass, including new parse/reject cases for `rate_correction`.

**Notes:**

Keep struct literals elsewhere compiling (tests, examples construct `SsaConfig { .. }`): add the field there or use `..Default::default()`.
