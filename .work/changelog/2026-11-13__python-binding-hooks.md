---
bump: minor
floor: minor
---
- Library entry points for the radiopyo Python binding: `Config::from_table` (a configuration from an already parsed TOML table), `GValuesBoscolo::from_csv_str` and `UnitCell::from_config_with_gvalues` (G-value data given by the caller, no file read).
