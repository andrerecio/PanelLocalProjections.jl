# Audit: the Almuzara–Sancibrián code

Almuzara & Sancibrián distribute `panel_LP` in MATLAB (`matlab_code/panel_LP.m`) and R
(`R_code/panel_LP.R`), at
[github.com/TinchoAlmuzara/PanelLocalProjections](https://github.com/TinchoAlmuzara/PanelLocalProjections)
(header "Version: 2024 Jun 10"; the audited copy equals the repository's `main` at commit
`4fd72b9`, 2026-08-16). The repository has **no licence**, so this package implements
the method from the paper's formulas only. The R function is *run* on this package's
synthetic panels, by `crosscheck/R/p4_as_panelLP.R`, and its outputs are committed as
test fixtures (`test/fixtures/r/p4_as_panelLP.csv`).

Statements marked *verified* were checked by running the R code on synthetic inputs;
the others come from reading the code.

## What the code computes

At each horizon ``h = 0, \dots, H``, `panel_LP` regresses ``y_{i,t+h}`` (or
``\sum_{j=0}^{h} y_{i,t+j}`` with `cumulative = TRUE`) on the interactions ``s_{it}
X_t``, ``p = \min(h, p_{\max})`` lags of ``y`` and of every ``s X`` column, and the
controls ``W``, after absorbing the fixed-effect columns passed in `FE` by alternating
projections (R: `fixest::demean`, tolerance 1e-8). The variance is clustered by period
with no small-sample factor and normal critical values; with `small_sample = TRUE`, HC2
on the synthetic regressors and Imbens–Kolesár degrees of freedom.

*Verified:* with `LagAugmented(CappedLags(p_max))`, `singletons = :keep` and
`TLAHR(weights = :raw)`, this package reproduces the code's estimates, standard errors
and degrees of freedom on 672 cross-check rows (five specifications × ``p_{\max} \in
\{0, 2\}`` × level/cumulative × with/without the refinement × a balanced and an
unbalanced panel): within 2e-12 on the balanced panel (estimates; standard errors 1e-13,
degrees of freedom 3e-13, relative) and 2e-9 on the unbalanced one, except for the bug B2
below.

## Differences from the paper

| | the code | the paper | `panellp` |
|---|---|---|---|
| lags | ``\min(h, p_{\max})``, ``p_{\max} = \lceil (T - H)^{1/3} \rceil``: none at ``h = 0``, no information criterion | ``\max\{1, \min(h, \hat p^{AIC})\}`` (eq. 15) | `ICLags()` by default; `CappedLags` reproduces the code |
| refinement | off by default (`small_sample = FALSE`) | recommended; "included by default in our own software" (supplement, fn 9) | `TLAHR()` refines by default |
| synthetic weights | the raw exposure ``s_{it}``, per period; periods unweighted | the demeaned exposure ``\hat s_i``, normalised by ``\sum_i \hat s_i^2`` (B.2: unbalanced panels are not covered) | residualised by default, periods weighted by ``\omega_t = \sum_i \hat s_{it}^2`` (fn 24; CR2 by period); `weights = :raw` reproduces the code |
| singletons | kept | — | dropped by default; `singletons = :keep` reproduces the code |
| without fixed effects | no intercept (a regression through the origin) | — | an intercept |

The two weightings coincide on balanced panels with time effects and a time-invariant
exposure (see [t-LAHR](../methods/tlahr.md)).

## Bugs

- **B1. Several aggregate shocks with the refinement crash.** The IK loop indexes the
  exposure matrix by the coefficient index, which runs past its columns when ``X`` has
  more than one column. *Verified:* with two columns in ``X`` and `small_sample = TRUE`,
  the R function stops with "subscript out of bounds"; without the refinement it runs.
- **B2. The degrees of freedom use the first coefficient's contrast for every
  coefficient** (`XX0[, 1]` in R, `XX0(:,1)` in MATLAB). *Verified:* with two exposures
  of one shock (``s`` and ``s^v``), the second coefficient's degrees of freedom equal the
  Imbens–Kolesár ν of its own synthetic regression computed with the *first*
  coefficient's contrast, to 6e-9. With its own contrast they differ by 0.06%–12% on our
  panels. `TLAHR()` uses each coefficient's contrast. A related design issue: with several
  exposures of the same aggregate series, the raw-exposure weights make the synthetic
  regressors nearly collinear. On our balanced panel, without lags (``h = 0``), ν is about
  4.5 for the first exposure, against 11.8 with weights residualised on the other exposure
  and 12.7 with that exposure alone. With two lags and a time-varying second exposure, the
  residualised weights do not remove the collinearity of the lags, and from ``h = 3`` ν is
  about 4.9 with them against 6.2–6.8 with the raw ones (see [t-LAHR](../methods/tlahr.md),
  The weights).
- **B3. Character fixed-effect columns silently change the estimates (R).** `cbind` of
  numeric and character columns coerces everything to character, `is.na("NaN")` is
  false, so rows with missing leads or lags stay in the sample, and their missing values
  are later set to zero. *Verified:* on the unbalanced panel, `FE = cbind(id, t)` as
  character gives 0.688 at ``h = 1`` against 0.600 with numeric columns (the ``h = 0``
  estimate, which needs no lead or lag, is unchanged).
- **B4. A `Date` time index breaks the function (R).** `as.matrix` strips the `Date`
  class before the date branch tests for it, so dates become day counts, and the grid step
  (the smallest gap, e.g. 90 days for quarters) leaves most periods off the grid.
  *Verified* on the unbalanced panel with quarterly dates: with `H = 0` the function
  returns 0.312 against 1.336 with integer periods, without a warning; with `H = 1` R
  crashes (a segmentation fault; the review traced it to `fixest::demean` on the empty
  sample of ``h = 1``).
- **B5.** The test for non-integer normalised periods is one-sided (`t − ⌊t⌋ > 0.01`),
  so a value just below an integer (e.g. 2.9999999999) is dropped, and the step is
  global: one off-grid date changes the grid for every observation.
- **B6. Documentation.** The README names the files `LP_panel.m`/`LP_panel.R` (they are
  `panel_LP.*`); the R usage example claims the same data-generating process as the
  MATLAB one, but draws the factor-loading noise from N(0, 1) instead of U(0, 1), with a
  different seed.

## Edge cases

- With ``s \equiv 1`` and time effects, ``X_t`` is absorbed, and the code returns a
  number. It warns (`SINGULARITY SAFEGUARD`) only when lags or controls are present, so at
  ``h = 0``, or with ``p_{\max} = 0``, the absorbed shock's estimate comes without a
  warning. `panellp` refuses an absorbed shock.
- Without fixed effects, no intercept is estimated.
- A period whose exposures are all zero (e.g. an empty group of a binary exposure) makes
  the refinement NaN for the whole horizon; `TLAHR()` gives that period zero weight.
- Duplicate (unit, period) rows are not detected; `panellp` refuses them.
