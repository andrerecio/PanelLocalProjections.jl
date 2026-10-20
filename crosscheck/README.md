# crosscheck/

Scripts that compute the reference numbers the Julia tests compare against: Stata through
PyStata, and R. The scripts are committed. The local toolchain and the run outputs
(`.venv/`, `crosscheck/ado/`, `crosscheck/out/`) are not.

## Setup (once per machine)

```bash
uv venv .venv --python 3.13
uv pip install --python .venv/bin/python -r crosscheck/requirements.txt
.venv/bin/python crosscheck/run_do.py crosscheck/stata/setup_stata.do --log crosscheck/out/setup_stata.log
.venv/bin/python crosscheck/ado_manifest.py --check     # installed files match the committed checksums
.venv/bin/python crosscheck/run_do.py crosscheck/stata/smoke.do
Rscript crosscheck/R/smoke.R
```

`setup_stata.do` installs every package into the project-local tree `crosscheck/ado/`.
`common.do` points `sysdir PLUS/PERSONAL` there, so the user's own Stata setup is never
touched. Versions and pins are in `VERSIONS.md`.

## Running a cross-check

```bash
.venv/bin/python crosscheck/run_do.py crosscheck/stata/<script>.do --log crosscheck/out/<script>.log
Rscript crosscheck/R/<script>.R
```

Each do-file starts with `do "$PLP_ROOT/crosscheck/stata/common.do"`. `run_do.py` sets
`PLP_ROOT` to the repository root and runs the do-file from its own directory, in a fresh
Stata session each time (~1–2 s). The exit code is the Stata return code.

## Conventions every script follows

- `set type double` (in `common.do`) and `import delimited …, asdouble`. Stata's float
  storage otherwise costs about 1e-7 relative error.
- Inputs are the committed synthetic panels in `test/fixtures/synthetic/`, or local-only
  data under `docs/references/` whose outputs go to `test/fixtures/local/` (gitignored).
- Outputs are written at full precision: Stata `%21x` or `%21.0g`; R `format(x, digits = 17)`
  or `data.table::fwrite`. Committed outputs go to `test/fixtures/{stata,r}/` and are
  listed with their provenance in `test/fixtures/README.md`.
- Re-runs are compared numerically (≤ 1e-12), ignoring header lines.
- fixest ≥ 0.13: `ssc(K.adj = , G.adj = )`, explicit `fixef.rm`, explicit `vcov =`.
- MSS `LP_panel` needs a base `data.frame` (`as.data.frame()` first). It is reliable only in
  its `H = 0` mode with pre-built leads and lags.

## Command ↔ Julia call map

Filled in phase by phase. Each row gives the exact `panellp`/`vcov` options that reproduce
a reference command: `singletons`, `timelags`, `ssc`, `psd`, `absorb`, `dist`.

Every row is a committed test (`test/test_crosscheck_p1.jl`), at horizons 0–6 on both
synthetic panels unless noted. The test requires each reference spec to have its exact
number of rows and a Julia mapping. It compares coefficients, standard errors,
observations, critical degrees of freedom, numbers of clusters or periods, and reghdfe's
`e(df_a)` and `e(df_a_initial)`. The formula is `leads(y) ~ x + lags(x, 1) + lags(y, 1) + w`
plus the fixed effects listed. "drop"/"keep" is the `singletons` keyword. Agreement is
the largest relative difference over all horizons, terms and panels.

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `xtlp F#.y …, fe m(fe) vce(cluster id)` | `fe(id)`, drop | `Clustered(:unit)` | b 5e-13, se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster id, ase)` | `fe(id)`, drop | `Clustered(:unit; ssc = :none)` (xtlp keeps t(G−1)) | se 6e-15 |
| `xtlp …, fe m(fe) vce(robust)` | `fe(id)`, drop | `Robust()` | se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster t)` | `fe(id)`, drop | `Clustered(:time)` | se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster id t [, ase])` | `fe(id)`, drop | `Clustered(:unit, :time [; ssc = :none])` | se 6e-15 |
| `xtlp …, fe m(fe) vce(dkraay [lag(3)])` | `fe(id)`, drop | `DriscollKraaySE(:nw94 \| 3; timelags = :position)` | se 6e-15 |
| `xtlp …, vce(dkraay lag(3), nodfadj \| ase)` | `fe(id)`, drop | `DriscollKraaySE(3; ssc = :xtscc \| :none, timelags = :position)` | se 5e-15 |
| `xtlp …, keepsingletons` | keep | as above | se 5e-15 |
| `xtlp …, tfe m(fe) vce(…)`, the same `vce()` options and `keepsingletons` | `fe(id) + fe(t)`, drop / keep | as above | b 7e-14, se 2e-11 (xtlp stops on a 1e-9 absolute step) |
| `xtlp y …, fe m(fe) hor(0 6) ytransf(level \| diff \| cmltdiff \| cmltsum)` | LHS `leads(y)`, `leads(firstdiff(y))`, `ldiff(y)`, `cumul(y)` | `Clustered(:unit)` | 1e-13 |
| `xtlp F#.y x w, fe \| tfe m(fe) vce(cluster id) [keepsingletons]` (no lags: units with one or two rows enter, and singleton drops cascade) | `leads(y) ~ x + w + fe(id) [+ fe(t)]` | `Clustered(:unit)` | se 3e-12 |
| `reghdfe …, absorb(id t) vce(cluster id t) tol(1e-14)` | `fe(id) + fe(t)`, drop | `Clustered(:unit, :time)` | se 5e-11 |
| `reghdfe …, absorb(id) vce(cluster reg)` | `fe(id)`, `clusters = (:reg,)` | `Clustered(:reg)` | se 9e-15 |
| `reghdfe …, absorb(id t grp#t) vce(cluster id)` | `fe(id) + fe(t) + fe(grp) & fe(t)` | `Clustered(:unit)` | se 2e-10 |
| `reghdfe …, absorb(id [t]) vce(dkraay 4)` | `fe(id) [+ fe(t)]` | `DriscollKraaySE(3; timelags = :position)` | se 5e-12 |
| `reghdfe …, absorb(id t) vce(robust)` | `fe(id) + fe(t)` | `Robust()` | se 2e-11 |
| `reghdfe y x w, absorb(id t) vce(cluster id)` (no lags) | `fe(id) + fe(t)`, drop | `Clustered(:unit)` | se 2e-11 |
| `reghdfe …, absorb(id) vce(cluster reg t)`: not PSD, the repair fires | `fe(id)` | `Clustered(:reg, :time)` | se 5e-15 |
| `reghdfe …, noabsorb vce(cluster reg t)`: pooled, `_cons` reported, the repair fires (slope block clipped on its own) | no fixed effects (intercept) | `Clustered(:reg, :time)` | se 1e-14 where reghdfe posts a V (*) |
| `xtreg …, fe vce(cluster id)`, also with `i.t` | `fe(id) [+ fe(t)]`, keep | `Clustered(:unit)` | se 2e-14 |
| `regress …, vce(cluster id)` | no fixed effects (intercept) | `Clustered(:unit)` | se 8e-15 |
| `xtscc …, fe lag(3)` (fx_balanced: xtscc refuses a period absent from the data for every unit) | `fe(id)`, keep | `DriscollKraaySE(3; ssc = :xtscc, timelags = :position)` | se 6e-15 |
| Derived panel `_t29` (t ≤ 29): `xtlp …, fe \| tfe m(fe) vce(dkraay)`, the default lags change with each horizon's T (3, then 2) | `fe(id) [+ fe(t)]`, drop | `DriscollKraaySE(:nw94; timelags = :position)` | se 2e-11 |
| Derived panel `fx_balanced_gap` (w missing at t = 20 for every unit, unit 30 at t = 10, 11 only): `xtscc …, fe lag(3)`; `xtlp …, fe m(fe) vce(dkraay lag(3))` | keep (xtscc), drop (xtlp) | `DriscollKraaySE(3; …, timelags = :position)` (calendar lags differ here) | se 6e-15 |
| Derived panel `fx_balanced_2blocks` (two connected components of unit and time effects): `reghdfe …, absorb(id t) vce(cluster id \| robust)` | `fe(id) + fe(t)`, drop | `Clustered(:unit)`, `Robust()` | se 4e-15; `e(df_a)` exact |
| R `feols(f(y,h) ~ … \| id [+ t], fixef.rm = "singleton" \| "none")` with `ssc(K.adj = FALSE, G.adj = FALSE)`: cluster, hetero, twoway | drop / keep | `Clustered(...; ssc = :none)`, `Robust(; ssc = :none)`, `Clustered(:unit, :time; ssc = :none, psd = :raw)` | 7e-14 |
| R `feols(...)`, default `ssc`, `cluster = ~id` | drop / keep | `Clustered(:unit)` | 2e-14 |
| R `feols(… \| id)` on units 1–4 (`fx_balanced_id4`), `cluster = ~id + t`: not PSD; default `vcov_fix = TRUE` / `vcov_fix = FALSE` | drop | `Clustered(:unit, :time; ssc = :none, psd = :raw \| false)` | 1e-10 (NaN where fixest's variance is negative) |

(*) At h = 0, 1 on fx_balanced Stata refuses reghdfe's repaired matrix when posting it
("variance matrix is nonsymmetric or highly singular") and reports zero standard errors.
The test then only requires that our repair fired. At h = 2 the repair fires and reghdfe
posts a valid matrix, which we match.

`xtlp …, tfe` counts the absorbed degrees of freedom as N + T − 1, which assumes the unit
and time effects form one connected component. `panellp` counts the components, as
reghdfe does. On a sample with disconnected blocks the two differ, so the `xtlp tfe`
mapping holds for connected samples only.
| R MSS `LP_panel(…, method = "FE", te, two_way_cluster, robust, eigens)` (H = 0, pre-built leads) | keep; `te = T` ⇒ `fe(id) + fe(t)`, `absorb = OneSweep()` | `Clustered(:unit [, :time]; ssc = :none \| :mss, psd = :none \| :raw)` | 7e-15 |
