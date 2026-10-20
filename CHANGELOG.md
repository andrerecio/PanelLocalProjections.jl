# Changelog

## Unreleased (0.1.0-DEV)

### Phase 3: cumulative-difference OLS (Ugarte-Ruiz 2026)
- `CDOLS(; validate = true)`: pooled OLS of the long difference y_{t+h} − y_{t−1} on the
  shock and its lag. The specification rules are checked structurally, whatever the
  spelling of the formula:
  - a long-difference response (`cumul(firstdiff(y))` is `ldiff(y)`);
  - every shock dated right after the base period, with its first lag among the
    regressors;
  - outcomes only lagged beyond the response's start, and differenced when the response
    differences them;
  - no fixed effect nested in the units, in the data or in any horizon's sample, with a
    warning for groups with a single unit.

  `validate = false` admits pre-built outcome columns.
- Cross-checked against `regress`, `xtscc` and `locproj … m(reg) tr(cmlt) sl(1) [yl(1)]`
  in Stata and fixest in R, on the panels drawn by `xtarsim` (the paper's Monte Carlo
  generator), the native design and an unbalanced panel. Coefficients agree within
  1e-12, standard errors within 5e-14, and observation, degrees-of-freedom and cluster
  counts exactly. The paper's comparators reproduce too: `xtreg, fe` (CMFE1) and
  `xtlp m(spj)` (the jackknife of its applications).
- Tests:
  - the probability limits of β_h, of the γ_h nuisance coefficient and of the response
    when s_{t−1} is omitted;
  - the different estimand of the levels LP without s_{t−1} (the paper's FE and
    jackknife comparators): +178% of β at h = 0 under joint persistence;
  - a negative test: unit-specific shock means make CD-OLS inconsistent;
  - that the native design, at the paper's signal-to-noise ratio, reproduces the
    `xtarsim` panels.
- A Monte Carlo of robust versus clustered standard errors in the paper's designs. With
  ρ_s = 0 and ρ_y = 0.8, robust ones overstate the sampling sd by 19–51% at h = 2–9 (the
  paper's over-coverage). With ρ_s = ρ_y = 0.8 they understate it by 7–12% at h = 1–4.
  Unit-clustered ones match it within 4%.
- The `xtarsim` panels are committed (kiss32 draws), with a test of their documented
  moments.
- The Cerra–Saxena application (local data) reproduces the paper's Fig. 9 CD-OLS path
  within 0.06. The Mian–Sufi–Verner application is not reproduced: the specification is
  not stated in the paper, and even its fixed-effects and jackknife lines differ from the
  replication file.
- Docs: `methods/cdols.md`, `replication/ugarte2026.md`, `audit/ugarte.md`.

### Phase 2: split-panel jackknife (Mei, Sheng & Shi 2026)
- `SplitPanelJackknife(; split = :ceil | :floor)`: per-unit halves of each horizon's
  sample, every fixed effect re-estimated inside each half, and the variance from the
  full-sample bread, the SPJ residuals and the scores d = 2X̃ − X̃_half.
- Cross-checked against:
  - `xtlp m(spj)`: 22 `fe`/`tfe` × `vce()` specs, including lag-free specs with
    single-row halves;
  - `xtspj` (Sun & Dhaene 2019);
  - the authors' `LP_panel`;
  - a hand-split identity and an explicit-dummy variance on an unbalanced panel, under
    both split rules.
- The authors' four applications are reproduced exactly: Table 1 to 3 decimals, and
  every horizon to about 1e-12. `xtlp tfe` confirms the exact two-way variant on the
  Romer–Romer data.
- Documented what the published numbers depend on:
  - The Romer–Romer unemployment path uses a one-sweep two-way demeaning, exact only on
    balanced samples. Exact two-way fixed effects lower the peak SPJ response by about 3%
    and its standard errors by 12–18%.
  - The authors' ⌈n/2⌉ split differs from the paper's ⌊n/2⌋. The paper's split moves
    the Romer–Romer GDP peak from −6.285 to −6.606.
  - The eigenvalue repair raises the Baron–Verner–Xiong standard errors by up to 13%.
  - `OneSweep()` and `split = :ceil` reproduce the published numbers.
- A regressor absorbed inside one half makes the horizon `:spj_half_singular` (NaN)
  instead of producing an arbitrary estimate.
- Monte Carlo test of the paper's toy design: FE Nickell bias, SPJ correction, coverage,
  and the calibration of the SPJ standard error.

### Phase 1: fixed-effects panel LP and inference
- `panellp(f, data; id, time, horizon, …)` with `WithinFE()`. Leads, lags and
  differences go through a (unit, period) lookup: they respect gaps and never cross units,
  and nothing is dropped before shifting. Integer or `Date` periods (with `frequency`).
  Per-horizon or common samples, a `subset` applied after the shifts, and iterative
  singleton removal.
- Formula helpers `leads`, `ldiff`, `anchor`, `cumul`, `lags`, `firstdiff`, `fe`
  (including `fe(a) & fe(b)`) and shock × exposure interactions; `@formula` re-exported.
- In-house fixed-effect absorption: exact one-way demeaning, alternating projections with
  an Irons–Tuck step, `OneSweep()` for the Mei–Sheng–Shi code. Absorbed degrees of freedom
  follow reghdfe's pairwise rule.
- OLS by column-equilibrated pivoted QR. Regressors absorbed by the fixed effects are
  omitted; an absorbed or collinear shock is an error.
- Post-estimation variances from stored score sums: `Robust`, `Clustered` (one- or
  two-way; `ssc = :none | :stata | :mss`; PSD repair `:stata | :raw`) and
  `DriscollKraaySE` (bandwidths `:nw94`, `:as2026`, `:mss2026` or fixed; calendar or
  position lags; `ssc = :none | :stata | :xtscc`). Critical values follow each convention.
- `summarize` → `IRFSummary` (a Tables.jl table), `coefpath`, `nobspath`, `esample`; a
  Plots recipe and a Makie extension (`irfplot`, `irfplot!`).
- Cross-checked against xtlp, reghdfe, xtreg, regress and xtscc (Stata), fixest and the
  Mei–Sheng–Shi `LP_panel` (R), at horizons 0–6 on the synthetic panels and on derived
  panels built for edge cases. These cover the default Driscoll–Kraay bandwidth per
  horizon, an interior period absent from every sample, singleton cascades, disconnected
  fixed-effect blocks, and two-way variances whose PSD repair fires. Coefficients agree
  within 5e-13 and standard errors within 2e-10. Observation and cluster counts, critical
  df and absorbed df agree exactly.
- Documentation: getting started, data, formula, inference conventions, absorption, API.

Fixes from the adversarial review of Phase 1:
- **Absorption:** the stopping rule of `AlternatingProjections` is relative to the column
  after its first sweep, so it no longer depends on the units or the level of the data.
  Small-scale columns were stopped early, and an aggregate shock at a small scale could
  escape the absorbed-shock error.
- **Dates:** `fe(date)` and clusters on a `Date` time column group by period, not by day.
- **Degenerate horizons:** a horizon whose sample leaves the shock absorbed or collinear
  is `:singular` (NaN, reported) instead of aborting every horizon. It is still an error
  when this happens at every horizon.
- **Undefined variances:** a single cluster or period, or no residual degrees of freedom,
  gives NaN with a warning instead of a LAPACK error. A negative variance gives a NaN
  standard error instead of a `DomainError`.
- **PSD repair (`psd = :stata`):** with an intercept, the slopes are clipped on their
  own and the intercept enters the metric with scale 1, as in reghdfe. fixest's default
  two-way repair is `psd = :raw`.
- **Missing values:** a function of a missing value is missing (`(w > 0)` was 0). `NaN`
  in a fixed-effect or cluster column is a missing value.
- **Warnings:** each call reports omitted regressors, empty or singular horizons and
  non-convergence. They are no longer suppressed after the first call of a session, and
  `show(r)` lists them.
- **Input checks:**
  - `ssc = :mss` only for `Clustered(:unit)` and `Clustered(:unit, :time)`;
  - two names for one cluster dimension, data columns that shadow `:unit`/`:time`, and
    unsupported time columns (which recursed endlessly) are errors;
  - `summarize(r, cv)` checks that `cv` was computed from `r`.
- **Inputs accepted:**
  - horizons as consecutive integer vectors;
  - `frequency` and `clusters` as strings;
  - numeric columns of abstract element type;
  - nested Bool functions (`ifelse`, `!`, `|`) in formulas;
  - formulas built in code from the helpers (they are StatsModels terms).

### Phase 0: bootstrap
- Package skeleton: `Project.toml` with Julia compat floor 1.10 (LTS) and versioned
  manifests (`Manifest-v1.12.toml`, `Manifest-v1.10.toml`, both gitignored).
- Test suite on TestItemRunner with a standalone `test/Project.toml` and Aqua (all checks
  enabled).
- Committed synthetic fixture panels (`test/fixtures/synthetic/`), generated by
  `dev/make_fixtures.jl` from the deterministic DGPs in `test/support/PLPDGPs.jl`.
- `dev/format.jl --check`, Documenter skeleton, GitHub workflows (CI on `lts` and `1`,
  check-only formatting, documentation, CompatHelper, TagBot, dependabot).
- Cross-check infrastructure in `crosscheck/`: PyStata runner, project-local Stata ado
  directory, pinned installs, smoke checks.
