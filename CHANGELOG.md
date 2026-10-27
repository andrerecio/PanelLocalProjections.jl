# Changelog

## 0.2.0 (2026-10-26): panel LP-IV

### Phase 6c: first-stage diagnostics, public example, docs
- `firststage(r; vcov = TLAHR(refine = false))` reports, per horizon, with any variance
  specification:
  - the Kleibergen–Paap (2006) rk Wald F, computed from the SVD of the normalised first
    stage, for one or several endogenous regressors;
  - the Montiel Olea–Pflueger (2013) effective F, for one endogenous regressor;
  - the first-stage coefficients and standard errors.

  The small-sample factor uses the first stage's degrees of freedom. The result is a
  Tables.jl `FirstStageSummary`. It reproduces ivreghdfe's `widstat` (just- and
  over-identified, one and three endogenous regressors, robust, one- and two-way clusters,
  Driscoll–Kraay) to 7e-11, and weakivtest's `F_eff` to 4e-14.
- New pinned Stata packages for the cross-checks: avar and distinct (needed by weakivtest).
- Review fixes:
  - the KP statistic uses the trailing singular vectors directly, an equivalent form of
    KP's eq. (12) that needs no inversion, so tied singular values (group dummies × an
    aggregate) no longer throw;
  - an indefinite two-way variance gives NaN with a warning, never a negative statistic.
    `psd` does not apply, as in ivreg2;
  - the first-stage scores are filled into one preallocated matrix.

  New Stata cases: m = 2 with L = 3, and a diagonal Π̂.
- Public-data example (`examples/us_states/`, docs "Examples"). A US state panel:
  - personal income (BEA);
  - exposure: the manufacturing share;
  - endogenous: the change in the 1-year rate (FRED GS1) × exposure;
  - instrument: the Jarociński–Karadi surprises × exposure.

  It shows that the first stage is weak clustered by period (KP F ≈ 2) and strong robust
  (≈ 35) or clustered by state (≈ 2·10⁵). It also shows the reduced form with refined
  t-LAHR. The scripts download and build the data into the gitignored `docs/data/`.
- Docs: `methods/lpiv.md`; the firm-level page maps CFFS23 and JL24 IV to `panellp` calls,
  cites the OW20 Econometrica replication package and marks the three studies as
  requiring WRDS.
- `.gitignore`: `docs/data/` (third-party replication packages and data built from them).
- Version 0.2.0.

### Phase 6b: lag augmentation and t-LAHR for panel LP-IV
- `LagAugmented` with an IV block (Almuzara & Sancibrián 2026, §3.4): at horizon h, p(h)
  lags of the outcome are exogenous, and p(h) lags of each endogenous shock are endogenous,
  instrumented by the same lags of every excluded instrument (`r.instlag`). Each horizon
  uses the instruments up to p(h).
- `ICLags` on IV fits selects the lags on the reduced form: the instruments take the place
  of the endogenous regressors, as in the OLS regression of the outcome on them. The
  instrument with the endogenous shock's exposure (`s & mp` for `s & ffr`) is the
  reduced-form shock and weights the criterion; the other instruments are plain regressors.
- Lag ℓ of the instruments enters only with lag ℓ of an endogenous shock, so the fit at
  horizon h depends only on p(h).
- `TLAHR(refine = false)` on these fits is A&S's LP-IV variance (supplement A.4).
  `TLAHR(refine = true)` remains an error on IV fits.
- Cross-checks of the lag-augmented regression: ivreghdfe (clustered by period, robust,
  `dkraay`) and fixest (clustered by period without correction, i.e. t-LAHR).
- Monte Carlo (`:mc`) of a panel LP-IV design with an endogenous aggregate and a noisy
  proxy (`PLPDGPs.lpiv_mc`):
  - the 90% t-LAHR intervals cover at 0.85–0.91, and unit clustering at 0.06–0.16;
  - IV removes OLS's +0.38 bias at h = 0.
- Adversarial review (code; econometrics with 10 mutations). Fixed:
  - the `ICLags` choice depended on the order of the instruments;
  - an absorbed extra instrument stopped `ICLags`;
  - instrument lags depended on pmax.

  Tests were added for instruments and outcomes that start later, for two lagged
  instruments, and for the fit notes; every mutation is now caught. Non-IV fits are
  bit-identical to 6a.

### Phase 6a: IV core
- Panel LP-IV (Almuzara & Sancibrián 2026, §3.4): a block `(endo ~ inst)` in the formula,
  e.g. `leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)`, makes each horizon's fit
  two-stage least squares on its fixed-effects design. The instruments are absorbed with
  the rest of the design, in the same buffer. The first stage is a Frisch–Waugh–Lovell
  projection on the exogenous regressors and the instruments, and the fitted values
  replace the endogenous columns. The fit stores the bread (X̂'X̂)⁻¹ and the scores X̂ ξ̂,
  where ξ̂ is the structural residual, so every variance applies unchanged:
  `Robust`, `Clustered`, `DriscollKraaySE`, `CosineHAR`, and `TLAHR(refine = false)`, which is
  the paper's LP-IV variance (supplement A.4). `TLAHR(refine = true)` on IV fits is an error
  for now.
- Several endogenous regressors are allowed (e.g. group dummies × an endogenous aggregate).
  The endogenous regressors are listed first and are the default `shock`.
- Instruments absorbed by the fixed effects, or collinear with the exogenous regressors or
  with earlier instruments, are left out of the first stage, with a note. A horizon with
  fewer instruments than endogenous regressors left, or whose instruments do not predict
  every endogenous regressor (a first-stage partial R² of round-off size, or collinear
  fitted values), is `:underidentified`. A model underidentified at every horizon is an
  error. IV fits store no synthetic series, so `synthetic()` on them is an error.
- Specification errors:
  - `fe()` or a constant inside the block;
  - two blocks;
  - an exogenous regressor listed as an instrument;
  - fewer instruments than endogenous regressors;
  - a nested block;
  - IV with `SplitPanelJackknife`, `CDOLS`, or `LagAugmented` (planned for 6b).
- Fixtures `fx_lpiv` and `fx_lpiv_unbal` (the paper's eq. 14 design).
- Cross-checks:
  - Stata: `ivreghdfe` (robust, one- and two-way clusters, `dkraay`), `xtivreg2, fe`, and
    several endogenous regressors;
  - R: fixest `feols` with an IV part.

  Coefficients agree to 5e-15 and standard errors to 2e-10 (unbalanced fixest: 5e-9).
  When a fixed effect is nested in the cluster, ivreghdfe's clustered small-sample factor
  counts two fewer absorbed degrees of freedom than reghdfe's (documented in
  `crosscheck/README.md`). New pinned Stata packages: ivreghdfe, ivreg2, ranktest,
  xtivreg2, weakivtest.

## 0.1.0 (2026-10-21)

### Phase 5: documentation, audit, release
- Docs: the cross-check map (generated from `crosscheck/README.md`), a table of the
  deviations from each reference implementation and the switches that reproduce it, audit
  pages for `xtlp`, tidyMacro's `fLPPanel` and Ramey's infrastructure code, and the
  firm-level applications (Ottonello–Winberry, Jeenas–Lagos) mapped to `panellp` calls.
- Memory: each fit absorbs the design in place in one buffer, of which the within
  transform, the kept columns and the scores are views; the raw columns are not kept (their
  statistics are taken first), and the responses are built per horizon. On the scale
  benchmark (1.1M rows, 21 horizons) the peak RSS falls from 1.6 to 1.3 GB (fixed effects),
  from 2.6 to 1.8 GB (`LagAugmented(ICLags())`) and from 3.2 to 2.2 GB (`CappedLags()`,
  up to 10 added lag columns). The estimates are unchanged.
- Every cross-check re-run from an empty `crosscheck/ado`: the installed Stata packages
  match the committed checksums, and every committed reference output is reproduced byte
  for byte.
- Version 0.1.0.

### Phase 4: t-LAHR inference (Almuzara & Sancibrián 2026)
- `LagAugmented(rule)`: at each horizon, p(h) lags of the outcome (y, Δy for `ldiff`, g
  for `cumul(g)`) and of every shock column are added to the fixed-effects regression,
  which needs unit effects (or finer ones nested in the units, e.g. firm × state).
  The rules:
  - `ICLags()`, the default: max{1, min(h, p̂)}, with p̂ chosen by AIC or BIC on the
    synthetic time series (the paper's eq. 15, fn 11 and supplement B.3–B.4). The choices
    the paper leaves open are pinned and documented.
  - `CappedLags(pmax)`: min(h, pmax), the authors' code.
  - `ConstantLags(p)`.
- `lagorder(r)` returns p(h); `r.lags`, a `LagSelection`, holds the criterion.
- `TLAHR(; refine = true, weights = :residual)`:
  - `refine = false` is CR0 clustered by period (eq. 11) with normal critical values.
  - `refine = true` adds HC2 on the synthetic regression and Student-t critical values
    with each coefficient's Imbens–Kolesár degrees of freedom (supplement B.2), from a
    thin SVD and a closed form without eigendecomposition or T × T matrix.
  - The synthetic weights are the exposure residualised per period (the paper's ŝ_i, and
    the exposures of the other regressors of the same aggregate series), and the
    synthetic regression is weighted by ω_t = Σ_i ŝ_it² (fn 24): with every regressor an
    aggregate series times the exposure, HC2 is then CR2 clustered by period, and ν its
    Bell–McCaffrey degrees of freedom. `weights = :raw` is the authors' code.
  - The exposure structure (an aggregate series times unit-level variables, lags
    included) is detected for every regressor, from the formula and the data. Each
    coefficient uses one synthetic regression: its shock's, its own, or the first
    shock's.
- `CosineHAR(B = :llsw)`: equal-weighted cosine variance of the period sums of the
  scores, as Stata's `harreg`.
- `synthetic(r, h)` returns the synthetic series of a shock.
- Cross-checked against:
  - the authors' `panel_LP.R`, run on our synthetic panels (672 rows): within 2e-12 on
    the balanced panel (estimates; standard errors 1e-13) and 2e-9 on the unbalanced
    one. The exception is their degrees-of-freedom bug with two exposures (e₁ contrast),
    which is verified and documented;
  - Stata: `reghdfe` for the lag-augmented regression (2e-14; the unbalanced panel's
    standard errors 5e-12), and `regress …, vce(hc2, dfadjust)` (HC2 and Bell–McCaffrey
    df, `e(adj_df)`) and `harreg` (EWC) on the synthetic regression (3e-14).
- Monte Carlo of the paper's stylized design (Table 1): coverage of unit, two-way,
  Driscoll–Kraay, plain and recommended t-LAHR intervals. At the paper's size (N = 1500,
  2000 replications, the `:mc_full` item) every row is within 0.04 of the published
  one. The recommended intervals cover 0.87–0.91 (published 0.87–0.90), and the
  footnote 10 calibration of σ_u² reproduces the unit-cluster row. A reduced version
  runs in the suite.
- Docs: `methods/tlahr.md` (with a sensitivity section), `audit/as_code.md`. The bugs of
  the authors' code are verified by running it on synthetic inputs.

Fixes from the adversarial review of Phase 4:
- **A crash in plain fixed-effects fits.** The synthetic weights solved a square system per
  period, which threw on valid data with two exposures of one shock (e.g. group dummies
  with a period of one unit), in every non-jackknife fit. The residualisation is now a
  rank-revealing projection, and a period whose residualised exposure is round-off has
  weight exactly zero. Round-off weights had made the refinement NaN at every horizon and
  corrupted the information criterion.
- **Label independence.** The weights are residualised on the exposures of every
  regressor of the same aggregate series, not only on the declared shocks, so the
  inference no longer depends on `shock`. Ownership is structural: `s2 & L1.X` uses the
  synthetic regression of the shock `s2 & X` (it had NaN degrees of freedom).
- **Unbalanced panels.** The synthetic regression is weighted by ω_t (above); the
  unweighted form gave small periods the leverage of full ones and overstated ν. The
  information criterion averages ξ̄_t² with the same weights, which the periods with few
  units dominated.
- **Lags.** A lag already among the regressors under another name (`s & L1.X` for the
  shock `s & L0.X`) is not added again: the duplicate had made horizons singular. A lag
  that the formula's columns span gives way to them in the collinearity screen, so the
  user's coefficients stay. The criterion counts only the lags a candidate adds, so lags of
  the outcome written in the formula change neither p̂ nor p(h).
- **Numerics and messages.** The degrees of freedom need memory linear in T (a T × T matrix
  took 2.9 GB at T = 8000). Specific messages for a cosine variance with B ≥ T, a
  `synthetic` call on a jackknife fit, a criterion sample emptied by `pmax`, and a shock
  absorbed by time effects under `ICLags`. A shock written `lags(s & X, 0:1)` is an
  exposure.
- **Tests.** An exact CR2 oracle, the zero-weight periods, the joint variance with its
  covariances, the default `pmax`, the criterion's variants on unbalanced data, and
  paired assertions in the Table 1 Monte Carlo that catch a missing half of the
  refinement.
- The scale benchmark (1.1M rows, 21 horizons) also runs `LagAugmented(ICLags())` with
  `TLAHR()`: 19 s and 0.7 s, 2.6 GB max RSS; and the authors' rule, `CappedLags()` (up to
  5 lags each of the outcome and of the shock): 25 s and 3.2 GB. The fixed-effects run:
  13.5 s, 1.6 GB. (Phase 5 lowers the memory.)

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
