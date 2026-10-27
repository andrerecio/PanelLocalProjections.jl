# Cross-check toolchain versions

Recorded on 2026-10-20 (Phase 0). Every committed reference output in `test/fixtures/` was
produced with the versions below, unless its header says otherwise. The checksums of the
installed Stata files are in `crosscheck/ado_manifest.sha256`. Verify them with
`.venv/bin/python crosscheck/ado_manifest.py --check`.

## Stata

StataNow 19.5 BE (`/Applications/StataNow`), driven by PyStata through `stata_setup`
0.1.3. The packages live in the project-local tree `crosscheck/ado/` (gitignored),
installed by `crosscheck/stata/setup_stata.do`.

| Package | Version (`*!` header) | Source | Pin / Distribution-Date |
|---|---|---|---|
| require | 1.4.0 (10apr2024) | GitHub sergiocorreia/stata-require | `7ee750b12f9c51ec8a879445cbf39357c308bdee` |
| ftools | 2.50.0 (09jan2026) | GitHub sergiocorreia/ftools | `7b3663e49ea5c5b81638c55be29edf416e68e8b7` |
| reghdfe | 6.14.1 (08Jul2026) | GitHub sergiocorreia/reghdfe | `29b2b203f534413369d675d0c0f674d61291f536` |
| xtlp | 1.1.1 | GitHub shenshuuu/panel-local-projection-stata | `cd3fa2468ada916d9c05323840b0a2abd48e0173` |
| xtscc | 1.4 (01dec2017) | SSC | 20180403 |
| locproj | version 10 (19 Feb 2026) | SSC | 20260219 |
| xtspj | (no header; compiled `lxtspj.mlib`) | SSC | 20190518 |
| xtarsim | 1.0.0 (31may05) | SSC | 20050620 |
| harreg | (no header) | SSC | 20260516 |
| ivreghdfe | 1.1.4 (29nov2025) | GitHub sergiocorreia/ivreghdfe | `bfb5577a6dbdfb029ab4ab6a7e93f7257a827b42` (added 2026-10-25, Phase 6) |
| ivreg2 | 4.1.12 (14aug2024) | SSC | 20240814 (Phase 6) |
| ranktest | 2.0.04 (21sept2020) | SSC | 20200929 (Phase 6) |
| xtivreg2 | 1.0.17 (19Feb2015) | SSC | 20200626 (Phase 6) |
| weakivtest | (no header) | SSC | 20260722 (Phase 6) |
| avar | 1.0.07 (28jul2015; a weakivtest dependency) | SSC | 20150730 (Phase 6c) |
| distinct | 1.2.1 (1mar2012; a weakivtest dependency) | SSC | 20120321 (Phase 6c) |

SSC packages cannot be pinned. A failed `ado_manifest.py --check` means an upstream
change: rerun the cross-checks, review the diff of the committed outputs, then
re-baseline with `ado_manifest.py`.

Built-ins used: `regress … vce(hc2, dfadjust)` (per-coefficient Bell–McCaffrey df in
`e(adj_df)`), `xtreg, fe`, `lpirf`.

## R

R 4.5.3 (2026-03-11): fixest 0.14.0, sandwich 3.1.1, data.table 1.18.2.1, pracma 2.4.6,
haven 2.5.5, lubridate 1.9.5, tidyverse 2.0.0. The reference code (local only) is
`docs/references/code/panel-lp-replication/applications/LP_panel_all.r` (MSS, MIT) and
`docs/references/code/PanelLocalProjections/R_code/panel_LP.R` (Almuzara–Sancibrián, no
licence; run only, never ported).

## Python

Python 3.13.12 in `.venv` (uv): stata_setup 0.1.3, numpy 2.5.4, pandas 3.0.6.

## Julia

Julia 1.12.6 (release) and 1.10.12 (LTS, the compat floor), with versioned manifests
`Manifest-v1.12.toml` and `Manifest-v1.10.toml`.
