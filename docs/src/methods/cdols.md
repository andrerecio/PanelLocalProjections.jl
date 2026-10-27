# Cumulative-difference OLS

Ugarte-Ruiz (2026) avoids the within transformation altogether. Taking the long difference
of the outcome removes the unit effect, so the local projection can be estimated by
pooled OLS:

```math
y_{i,t+h} - y_{i,t-1} = c_h + \beta_h s_{i,t} + \gamma_h s_{i,t-1} + \varepsilon_{i,t+h},
\qquad h = 0, \dots, H .
```

No unit effect is estimated, so the incidental-parameter (Nickell) bias of the
fixed-effects estimator does not arise. Because ``y_{i,t-1}`` is predetermined,
``\beta_h`` is the response of the *level* ``y_{i,t+h}``: the estimand of a fixed-effects
LP in levels that also controls for ``s_{i,t-1}``. A levels LP *without* the shock's lag
(the paper's fixed-effects and split-panel jackknife comparators, its eq. 24) targets
another quantity when both series are persistent: it converges to ``\beta_h + \beta
\rho_y^{h+1}\rho_s/(1 - \rho_y\rho_s)``, +178% of ``\beta`` at ``h = 0`` with
``\rho_y = \rho_s = 0.8`` (a test item checks both limits on a long panel).

```julia
r = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), data; id = :id, time = :t,
            horizon = 9, estimator = CDOLS())
summarize(r, Clustered(:unit))
```

## The lag of the shock is required

With a persistent shock (``\rho_s \neq 0``), ``s_{i,t}`` is correlated with the
unobserved past of ``y`` that remains in the long difference. Partialling out
``s_{i,t-1}`` leaves the shock innovation, which is orthogonal to it. In the AR(1) design
of the paper (``y_{it} = \rho_y y_{i,t-1} + \beta s_{it} + u_i + e_{it}``,
``s_{it} = \rho_s s_{i,t-1} + v_{it}``):

| quantity | probability limit |
|---|---|
| ``\hat\beta_h`` | ``\beta_h = \beta \sum_{j=0}^{h} \rho_y^{h-j} \rho_s^{j}`` |
| ``\hat\gamma_h`` (a nuisance coefficient) | ``(\rho_y^{h+1} - 1)\beta / (1 - \rho_y \rho_s)`` |
| ``\hat\beta_h`` without ``s_{i,t-1}`` | ``\beta_h - (1 - \rho_y^{h+1}) \beta \rho_s / (1 - \rho_y \rho_s)`` |

The test suite checks all three on large simulated panels. Omitting ``s_{i,t-1}`` biases
the response by −80% of ``\beta`` at every horizon when only the shock is persistent
(``\rho_s = 0.8``, ``\rho_y = 0``), and by −44% at ``h = 0`` to −198% at ``h = 9`` when
both are (``\rho_s = \rho_y = 0.8``; −47% to −309% at 0.9); not at all when
``\rho_s = 0``. `CDOLS()` therefore requires it.

## Specification rules

`CDOLS()` checks the formula (`CDOLS(validate = false)` skips the checks, e.g. for a
pre-built outcome column):

- the response is the long difference `ldiff(y)`, or `cumul(g)` of a growth rate `g`
  (``\sum_{k=0}^{h} g_{i,t+k} = y_{i,t+h} - y_{i,t-1}``); `cumul(firstdiff(y))` is
  `ldiff(y)`;
- every shock is dated ``t``, right after the base period ``t - 1`` (with a lagged
  response, `ldiff(lags(y, ℓ))`, it is dated ``t - \ell``): a shock dated ``t - 1``
  against the base ``y_{i,t-1}`` would estimate ``\beta_{h+1} - \beta_0``;
- the first lag of every shock is a regressor, however it is written (`lags(s, 1)`,
  `lags(s, 0:1)` with `shock = "L0.s"`, `lags(firstdiff(s), 1)`, a pre-built column, …);
- outcome-derived regressors are lagged beyond the response's start, and an outcome that
  the response differences enters only in differences, `lags(firstdiff(y), 1:p)` (lags of
  a growth rate `g` are allowed). A *level* ``y_{i,t-1}`` contains the unit effect
  ``u_i/(1 - \rho_y)``: it biases ``\hat\beta_h`` when the unit effects are correlated
  with the shock, as in the paper's design;
- no unit fixed effect, nor one nested in the units, in the data or in any horizon's
  sample (e.g. after `subset` or missing values): within estimation would bring the
  Nickell bias back. Time effects, `fe(t)`, are allowed. A grouping with few units per
  group brings back part of that bias (about ``1/g`` of it with ``g`` units per group),
  so a warning reports groups that contain a single unit.

The paper also gives every control the same lags as the shock (a contemporaneous
persistent control without its lag biases the response, like an omitted
``s_{i,t-1}``); that is not checked. `0 +` drops the intercept, as `noconstant` in
`locproj`.

## Inference

The paper reports "robust" standard errors without a formula. The long-difference error
is serially correlated, so heteroskedasticity-robust standard errors are not valid.
Simulating the paper's designs (``N = 100``, ``T = 50``, signal-to-noise ratio 20; 2000
panels):

- with ``\rho_s = 0``, ``\rho_y = 0.8`` (its Table 5), robust standard errors overstate the
  sampling standard deviation of ``\hat\beta_h`` by 19%, 51% and 34% at ``h = 2, 4, 9``,
  and robust 95% intervals cover 98–99.8% there: the paper's over-coverage;
- with ``\rho_s = \rho_y = 0.8`` (its main design), they understate it by 7–12% at
  ``h = 1``–``4`` (coverage about 93%), and overstate it by 12% at ``h = 9``.

Standard errors clustered by unit match the sampling standard deviation within 4% in
both designs, with about 95% coverage. Clustering by unit is valid for a fixed ``T``,
Driscoll–Kraay for a large ``T``. A test item checks these orderings with 200 panels.

| reproduces | `vcov` spec |
|---|---|
| `regress …, vce(robust)`; `locproj …, vce(robust)` (the paper, probably) | `Robust()` |
| `regress …, vce(cluster id)`; `locproj …, vce(cluster id)`; fixest's default | `Clustered(:unit)` |
| `xtscc …, lag(3)` | `DriscollKraaySE(3; ssc = :xtscc, timelags = :position)` |

## When it fails

The probability limits above assume what the paper's design builds in:

- a common mean and a common ``\rho_s`` of the shock across units. With unit-specific
  shock means, the pooled regression of ``s_{i,t}`` on ``s_{i,t-1}`` leaves the mean in
  the residual, and CD-OLS is inconsistent. The test suite demonstrates a large bias;
- ``|\rho_y| < 1`` with a stationary start. With a unit root or unit-specific trends, the
  long difference keeps a unit term ``(h+1)g_i``;
- no feedback from ``y`` to ``s``.

The paper's Proposition 1 shows consistency, not finite-sample unbiasedness (see the
[audit](../audit/ugarte.md)).

## Cross-checks

The tests compare `CDOLS()` with Stata (`regress`, `xtscc`, `locproj …, m(reg)
tr(cmlt) sl(1) [yl(1)]`) and fixest, at horizons 0–9 on five panels:

- the native design `fx_cdols` (``\rho_y = \rho_s = 0.8``, correlated effects);
- three panels drawn by `xtarsim` (Bruno 2005), the generator of the paper's Monte Carlo,
  with ``(\rho_s, \rho_y) \in \{(.8, 0), (0, .8), (.8, .8)\}``, ``N = 100``, ``T = 50``,
  signal-to-noise ratio 20 and the paper's correlated effects
  ``u_i = 10(1 - \rho_y)(1 + \bar s_i - \bar s)``. The panels are committed, because the
  draws come from Stata's kiss32 generator. A test checks their moments against
  `xtarsim`'s formulas, and another that the native design, at the same signal-to-noise
  ratio, reproduces them;
- `fx_unbalanced` (shock `x`): gaps, singleton units, and a period absent for every
  unit (so no `xtscc` there, which refuses it).

Coefficients agree within 1e-12, standard errors within 5e-14, and observation counts,
degrees of freedom and cluster counts exactly. The paper's comparators reproduce as
well: `xtreg, fe` on the long difference (the paper's CMFE1), and `xtlp …, m(spj)` on
the level response (the split-panel jackknife of the paper's applications; its Monte
Carlo used `locproj …, m(xtspj)`, which averages two splits when a unit has an odd number
of rows and reports a within-OLS variance with normal critical values).

| reference command | `panellp` call |
|---|---|
| `regress F#.y−L.y s L.s, vce(cluster id \| robust)` | `ldiff(y) ~ s + lags(s, 1)`, `CDOLS()`; `Clustered(:unit)` / `Robust()` |
| `… LD.y`, `… L2.s` | `+ lags(firstdiff(y), 1)`, `lags(s, 1:2)` |
| `… i.t` | `+ fe(t)`, `singletons = :keep` (`regress` keeps singleton periods) |
| `xtscc F#.y−L.y s L.s, lag(3)` | `DriscollKraaySE(3; ssc = :xtscc, timelags = :position)` |
| `locproj y s, m(reg) tr(cmlt) sl(1) [yl(1)] vce(…)` | as `regress` (`yl(1)` is `lags(firstdiff(y), 1)`) |
| fixest `feols(ld ~ s + l(s, 1))`, `cluster = ~id` / `"hetero"` | `Clustered(:unit)` / `Robust()` |
