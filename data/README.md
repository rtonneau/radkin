# Bundled data

Copied unchanged from radiopyo (https://github.com/rtonneau/radiopyo, commit e0c053d,
LGPL-3.0, author R. Tonneau):

| File | radiopyo source |
|---|---|
| `advanced_water_radiolysis.toml` | `radiopyo/data/advanced_water_radiolysis.toml` |
| `reactions_biology.toml` | `radiopyo/data/reactions_biology.toml` |
| `config_Labarbe.toml` | `radiopyo/data/config_Labarbe.toml` |
| `gvalues_from_Boscolo_2020.csv` | `radiopyo/data/G_VALUES/gvalues_from_Boscolo_2020.csv` (Boscolo et al. 2020, doi:10.3390/ijms21020424) |
| `CONFIG_PAPER_01.toml` | `examples/tonneau_2025_pmb/data/CONFIG_PAPER_01.toml` (Tonneau et al. 2025, PMB 70 235021) |

`tonneau_2025.toml` is radkin's own: it includes `CONFIG_PAPER_01.toml` and adds the inputs
radiopyo's `examples/tonneau_2025_pmb/paper_runs.py` sets in code for the original runs
(variable G-values, LET 0.56 keV/µm, a 30 Gy single pulse at 1e7 Gy/s).

Include entries starting with `radiopyo/` (or `radkin/`) in a config resolve to this folder.
