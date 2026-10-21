# t-LAHR inference

Almuzara & Sancibrián (2026) study panel local projections of micro outcomes on a *macro*
shock ``X_t`` interacted with a unit-level exposure ``s_i``:

```math
y_{i,t+h} = \beta(h)\, s_i X_t + \psi(h)' W_{it} + \mu_i(h) + \nu_t(h) + \xi_{it}(h),
```

with unit and time effects (time effects only when ``s_i`` varies; with ``s_i \equiv 1``
they would absorb ``X_t``). The estimand is the slope of the best linear predictor of the
unit responses on the exposure, ``\beta(h) = \operatorname{Cov}(s_i, \beta_{ih}) /
\operatorname{Var}(s_i)``. Aggregate shocks make the errors dependent across units, so the
effective sample size is ``T``, not ``NT``. Their recommendation, **t-LAHR**, is

- **t**ime-clustered: the variance is clustered by period (CR0, eq. 11);
- **l**ag-**a**ugmented: the regression includes ``p(h)`` lags of the outcome and of the
  shock, ``W_{it} = (s_i X_{t-1}, y_{i,t-1}, \dots, s_i X_{t-p}, y_{i,t-p}, \tilde W_{it})``.
  With ``p \ge h`` the scores ``X_t \xi_t`` are serially uncorrelated (their Lemma A.1), so
  no kernel or bandwidth is needed;
- **h**eteroskedasticity-**r**obust: the variance is the Eicker–Huber–White variance of
  the *synthetic* time-series regression (Remark 1), refined in small samples with HC2 and
  Student-t critical values with the Imbens–Kolesár degrees of freedom (supplement B.2).

It is not a HAR estimator: Driscoll–Kraay, equal-weighted cosine (EWC) and clustering by
unit or two-way are the paper's comparators.

```julia
r = panellp(@formula(leads(y) ~ s & X + fe(id) + fe(t)), data; id = :id, time = :t,
            horizon = 12, estimator = LagAugmented())       # ICLags(): eq. (15)
summarize(r, TLAHR(); level = 0.90)                         # HC2 + Imbens–Kolesár t(ν)
lagorder(r), r.lags                                         # p(h); the criterion
synthetic(r, 4)                                             # the synthetic series at h = 4
```

## Lag augmentation

[`LagAugmented`](@ref) adds, at each horizon, ``p(h)`` lags of the outcome and of every
shock column (named `L1.y`, `L1.(s & X)`, …). A lag already among the regressors is not
added again, whatever its name: columns are compared numerically, so with the shocks
`s & lags(X, 0:2)` the lag of `s & L0.X` is `s & L1.X`. A lag that the formula's columns
span is left out of the regression, without a warning: with `lags(y, 1:2)` in the formula,
the lag ``\Delta y_{t-1} = y_{t-1} - y_{t-2}`` that `ldiff(y)` adds is redundant, and the
formula's two coefficients stay. A lag of a shock lags the whole column, so a
time-varying exposure enters as ``s_{i,t-\ell} X_{t-\ell}``, as in the authors' code (the
paper writes ``s_i X_{t-\ell}``). The lagged outcome depends on the response:

| response | lagged outcome |
|---|---|
| `leads(y)`, `anchor(y, z)` | ``y`` |
| `ldiff(y)`, `leads(firstdiff(y))` | ``\Delta y`` |
| `cumul(g)` | ``g`` (the authors' `cumulative = true`) |

The rule ``p(h)``:

| rule | ``p(h)`` | source |
|---|---|---|
| [`ICLags`](@ref)`()` (default) | ``\max\{1, \min(h, \hat p)\}``, ``\hat p`` by AIC | the paper, eq. (15) and B.7 |
| `ICLags(criterion = :bic)` | the same with BIC | supplement B.3 |
| [`CappedLags`](@ref)`(pmax)` | ``\min(h, p_{\max})``: no lag at ``h = 0`` | the authors' `panel_LP.m`/`.R` |
| [`ConstantLags`](@ref)`(p)` | ``p`` | "plain" t-LAHR (``p = 1``, Table 1) |

The paper keeps at least one lag at every horizon, including ``h = 0``, as insurance
against near unit roots.

### The information criterion

The paper chooses ``\hat p`` on the synthetic time series (fn 11; supplement B.1): a panel
criterion with penalty in ``(NT)^{-1}`` over-fits, because redundant lags improve the fit
at the rate ``T^{-1}``. The paper leaves several details open; `ICLags()` pins them:

- the "panel VAR(p)" is the ``h = 0`` regression of the outcome on the formula's
  right-hand side (with the contemporaneous shock) and ``p`` of its own lags, with the
  formula's fixed effects: B.3 adds one parameter per lag, the lag of ``\bar Y``
  (`vars = :outcome_and_shock` also adds the shock's lags, with penalty
  ``p(1 + n_{\text{shock}})g``);
- each candidate ``p`` is fitted on its own sample, the rows where its lags are defined
  (B.4 sums over ``t = p+1, \dots, T``; `sample = :common` uses the ``p_{\max}``-lag
  sample for all);
- the residuals are averaged per period with the synthetic weights of the first shock
  (below), ``\bar\xi_t(p)``, and ``\hat\sigma^2(p) = \sum_t \omega_t \bar\xi_t(p)^2 /
  \sum_t \omega_t`` over the periods with positive weight ``\omega_t``: the mean square
  residual of the synthetic regression, which fn 24 weights by ``\omega_t``. On a balanced
  panel with a time-invariant exposure ``\omega_t`` is constant and ``\hat\sigma^2`` is the
  mean of ``\bar\xi_t^2``, as in B.4. Unweighted, the periods with few units, whose
  ``\bar\xi_t`` have variance ``\propto 1/\omega_t``, would dominate the criterion (with
  staggered entry it then often selects no lag, or ``p_{\max}``, for a strong AR(2));
- ``\mathrm{IC}(p) = \log\hat\sigma^2(p) + k(p)\,g``, with ``k(p)`` the lags candidate
  ``p`` adds to the regression, ``g = 2/T`` (AIC) or ``\log T / T`` (BIC) and ``T`` the
  periods of the ``p = 0`` regression. A lag already in the formula (or spanned by it)
  adds nothing: candidates that differ only by such lags are the same regression, and
  ``\hat p`` is the last of them, so writing `lags(y, 1)` in the formula changes neither
  ``\hat p`` nor ``p(h)``. Other ties go to the smaller ``p``;
- ``p_{\max}`` is not given in the paper. The default is ``\max\{1, \lceil (T -
  H)^{1/3} \rceil\}``, the authors' code's ``p_{\max}`` (``T`` the distinct periods in the
  data, ``H`` the largest horizon): with ``T = 124`` and ``H = 20`` it is 5; the paper's
  applications select up to 4.

The criterion needs the first shock to be an aggregate series times an exposure; for
other shocks use `ConstantLags` or `CappedLags`.

### Sensitivity

On the cross-checks' panels (`fx_tlahr`, `fx_tlahr_unbal`: a first-order design) every
option of the criterion (`criterion = :bic`, `sample = :common`,
`vars = :outcome_and_shock`, `pmax` 2 or 6) selects ``\hat p = 1``, as `ICLags()`; on the
test suite's AR(2) panel (`ar2panel()` in `test/test_tlahr.jl`: ``N = 200``, ``T = 60``,
autoregressive coefficients 1.2 and −0.5) every option selects ``\hat p = 2``. What moves
the inference is the rule. Standard error of the shock, with ``\nu`` in parentheses, at
``h = 0, 2, 4``:

| rule (``p(h)`` for ``h = 0, 1, 2, \dots``) | `fx_tlahr_unbal` | AR(2) panel |
|---|---|---|
| `ICLags()`, any option (1, 1, 1, … and 1, 1, 2, 2, …) | .110 (11.9), .207 (11.4), .255 (10.6) | .060 (15.9), .167 (13.4), .263 (13.2) |
| `CappedLags()` (0, 1, 2, 3, 4, 4, …) | .167 (12.6), .251 (9.3), .311 (5.7) | .229 (16.1), .167 (13.4), .278 (12.9) |
| `ConstantLags(1)` (1, 1, 1, …) | as `ICLags()` | .060 (15.9), .178 (15.3), .260 (15.0) |

Without a lag at ``h = 0`` (the authors' code) the standard error rises by half on
`fx_tlahr_unbal` and almost fourfold on the AR(2) panel: the lagged outcome absorbs the
persistence of ``y``. At long horizons the code adds ``p_{\max} = 4`` lags each of the
outcome and of the shock, eight synthetic regressors, and ``\nu`` falls to 5.7. One lag
where two are needed (`ConstantLags(1)` on the AR(2) panel) moves the standard error at
``h = 2`` by 7% and raises ``\nu``.

## The variance

### CR0 by period

`TLAHR(refine = false)` is the variance clustered by period with no small-sample factor
and normal critical values, eq. (11):

```math
\hat\sigma^2(h) = \frac{\sum_t \big(\sum_i \tilde x_{it}\hat\xi_{it}\big)^2}
                       {\big(\sum_{i,t} \tilde x_{it}^2\big)^2},
```

``\tilde x_{it}`` being the shock residualised on the other regressors and the fixed
effects. It equals `Clustered(:time; ssc = :none)`. With [`WithinFE`](@ref) and no lags it
is the paper's "t-HR" comparator.

### HC2 and Imbens–Kolesár degrees of freedom

`TLAHR()` transfers the small-sample refinements of heteroskedasticity-robust inference
to the synthetic regression (supplement B.2). For a shock with weights ``w_i``:

1. the synthetic regressors are the per-period weighted averages of the
   within-transformed regressors, ``\bar X_t = \sum_{i \in t} w_i \ddot X_{it} /
   \omega_t``, ``\omega_t = \sum_{i \in t} w_i^2``, one row per period of the sample;
2. the synthetic regression is weighted by ``\omega_t`` (fn 24), so its rows are
   ``\tilde X_t = \sqrt{\omega_t}\,\bar X_t``, and ``d_t = 1 - h_t``, with ``h_t`` the
   leverage of period ``t`` in ``\tilde X``;
3. the variance is ``B\,(\sum_t S_t S_t' / d_t)\,B``, with ``S_t`` the scores summed in
   period ``t`` and ``B`` the bread of the panel regression;
4. the critical value is ``t_{1-\alpha/2}(\nu_k)``, with the Bell–McCaffrey degrees of
   freedom of coefficient ``k`` under a homoskedastic reference (Imbens & Kolesár 2016):
   ``\nu_k = (\sum\lambda)^2 / \sum\lambda^2``, ``\lambda`` the eigenvalues of ``G'G``,
   ``G_{\cdot t} = \bar M_{\cdot t} a_t / \sqrt{d_t}``, ``a = \tilde X (\tilde X'\tilde
   X)^{-1} e_k``, ``\bar M = I - \tilde X (\tilde X'\tilde X)^{-1} \tilde X'``. Since
   ``G'G = D \bar M D`` with ``D = \operatorname{diag}(a / \sqrt d)``, ``\sum\lambda =
   \sum_t a_t^2`` and ``\sum\lambda^2 = \sum_{t,s} a_t^2 a_s^2 \bar M_{ts}^2 / (d_t
   d_s)``, computed without an eigendecomposition.

**The period weights.** B.2 normalises the averages by ``\sum_i \hat s_i^2``, which is
constant over time on a balanced panel, and does not say what to do when ``\omega_t``
varies (unbalanced panels, time-varying exposures). When every regressor is an aggregate
series times the exposure, ``\ddot X_{it} = w_i \bar X_t'``, the block of period ``t`` in the
panel's hat matrix has rank one, and its eigenvalue is the leverage of
``\sqrt{\omega_t}\,\bar X_t``: the HC2 variance of the weighted synthetic regression is then
CR2 clustered by period (Bell & McCaffrey 2002), and ``\nu_k`` its Bell–McCaffrey degrees
of freedom. A test checks both on an unbalanced panel, to 1e-12. The authors' code uses
the unweighted ``\bar X_t``, which gives a period with few units the leverage of a full
one; `TLAHR(weights = :raw)` reproduces it. On a balanced panel with a time-invariant
exposure ``\omega_t`` is constant and the two coincide.

**Stata.** ``\nu_k`` is, by construction, the degrees of freedom of Stata's `regress …,
vce(hc2, dfadjust)` (`e(adj_df)`) on the weighted synthetic regression. The standard
error is the panel's, ``B\,(\sum_t S_t S_t'/d_t)\,B``; it equals that regression's HC2
standard error only when the panel regression is itself a synthetic regression (every
regressor an aggregate series times the exposure, on a balanced panel), as in the
cross-checks below. With the lagged outcome that `LagAugmented` adds, a unit-level
regressor, they differ by a few percent (the paper's Remark 1 says "≈").

**Numerics.** ``\tilde X'\tilde X`` is never formed: ``\tilde X`` has a thin SVD, and
singular values below ``\max(T, K)\,\varepsilon\,\sigma_1`` are dropped (a comment in the
authors' code notes that ``\bar X'\bar X`` can be near-singular, ``\kappa \approx
10^{11}``, and that a looser pseudo-inverse then makes the degrees of freedom depend on the
truncation rank). With ``U`` its left singular vectors and ``b = a^2/d``, ``\sum\lambda^2 =
\sum_t b_t^2 d_t^2 + \lVert U' \operatorname{diag}(b)\, U \rVert_F^2 - \sum_t b_t^2 h_t^2``:
no ``T \times T`` matrix is formed, so the memory is linear in the number of periods
(daily panels of thousands of periods). A period with leverage one (as many periods as
synthetic regressors, without unit effects) makes the variance undefined: NaN, with a
warning.

**The absorbed intercept.** B.2 builds ``\bar X`` from the two-way within-transformed
regressors, without a constant. On a balanced panel with time-invariant weights and unit
effects, the columns of ``\bar X`` sum to zero over the periods, so the leverages add up
to at most ``\operatorname{rank}(\bar X) \le T - 1`` and every ``d_t \ge 1/T``: the
synthetic regression behaves as one with an intercept whose leverage is not counted.

### The weights

The paper's weights are the demeaned exposure, ``\hat s_i = s_i - \bar s``, period by
period when ``s_{it}`` varies (Remark 3). `TLAHR()` uses, for a regressor ``s_k X_t``, its
exposure residualised per period

- on an intercept, when a fixed effect nested in the periods is absorbed (time effects,
  or group × time effects, within whose groups it is demeaned), and
- on the exposures of the other regressors that multiply the same aggregate series (up
  to scale), shocks or not.

The second rule extends the paper, which covers a scalar exposure. With several exposures
of one shock (e.g. turnover in two liquidity groups, as in Jeenas & Lagos 2024), the raw
exposures make ``\bar X`` nearly singular. By Frisch–Waugh–Lovell, the coefficient on
``s_k X_t`` is identified by ``s_k`` net of the other exposures of ``X_t``. The panel's
residual is global, with a coefficient weighted by ``X_t^2``; the residual period by
period equals it on a balanced panel with time-invariant exposures, where the other
exposures' synthetic columns then vanish (up to rounding, which the rank tolerance
drops). Since the rule follows the regressors, not their labels, the inference on
``s_k X_t`` does not depend on which of them `shock` names. With time-varying exposures or
an unbalanced panel the other columns do not vanish: the contemporaneous ones survive
through the unit effects, and the lags of the other exposure that `LagAugmented` adds stay
nearly collinear with the shock's own lags. ``\bar X`` is then ill-conditioned and ``\nu`` low: with two lags and a
time-varying second exposure on the cross-checks' balanced panel, ``\nu`` of the first
exposure is about 4.9 from ``h = 3`` on, against 6.2–6.8 with the raw exposures and
8.1–9.6 with that exposure alone (at ``h = 0``, without lags: 11.8, 4.5 and 12.7). These
intervals err on the wide side. A period whose residualised exposure vanishes (no more
units than exposures, or one exposure common to all units under time effects) has weight
exactly zero, not a round-off weight.

`TLAHR(weights = :raw)` weights by the raw exposure ``s_{it}`` and leaves the periods
unweighted, as the authors' code. The two agree on balanced panels with time effects and a
time-invariant exposure (where the within-transformed regressors sum to zero in each period,
the weights differ by a constant and ``\omega_t`` is constant), and for ``s \equiv 1`` on
balanced panels. They differ on unbalanced panels, with time-varying exposures, and with
several exposures. On the cross-checks' unbalanced panel, the shock's standard errors
differ by up to 3% and its degrees of freedom by up to 16%; with a time-varying exposure,
by up to 4% and 29% (1% and 5% on the balanced panel).

### Several shocks and the other coefficients

Every coefficient takes its ``d_t`` and degrees of freedom from one synthetic regression:

- a shock, its own, as in the authors' code (which reports one variance per shock);
- a lag of a shock (`L1.(s & X)`, however it was added), and a regressor with the same
  exposure as a shock (`s & L1.X` for the shock `s & X`), the shock's;
- any other aggregate series times an exposure (e.g. a control `s2 & X`), its own;
- every other coefficient (outcome lags, unit-level controls), the first shock's.

Writing ``z^{(j)}_t = (B S_t)_j / \sqrt{d^{(c(j))}_t}`` for coefficient ``j`` with owner
``c(j)``, the variance is ``V_{jl} = \sum_t z^{(j)}_t z^{(l)}_t``. Each coefficient's own
variance is the paper's HC2 variance, and ``V`` is positive semi-definite. The degrees of
freedom of coefficient ``j`` use its own contrast ``e_j``; the authors' code uses ``e_1``
for every coefficient (see the [audit](../audit/as_code.md)).

`TLAHR()` needs every shock to be an aggregate series (equal across units in every
period), or an interaction of aggregate series with unit-level variables: `s & X`,
`X & s`, `log(X) & s & w`, `lags(s & X, 0:1)`, `X` alone (exposure 1). For any other
shock, `TLAHR()` is an error and `TLAHR(refine = false)` gives the CR0 variance. Neither is
defined for `SplitPanelJackknife` results.

### Unit-level controls

B.2 averages the whole regressor vector, so each regressor enters the synthetic regression
and takes up leverage there, unit-level controls included, although the panel identifies
their coefficients from cross-sectional variation. Irrelevant controls therefore widen
the refined intervals. On `fx_tlahr` with one lag, twenty columns of noise raise the
shock's refined standard error by about two thirds, at ``h = 0`` (0.096 to 0.158) and at
``h = 4`` (0.256 to 0.421), and lower ``\nu`` at ``h = 4`` from 10.4 to 7.5, while the
coefficient and the CR0 standard error do not move. Include the controls the specification needs, no more;
`TLAHR(refine = false)` does not depend on them.

## The synthetic series

[`synthetic`](@ref)`(r, h)` returns, per period of the horizon's sample, the weighted
averages of the within-transformed shock (`x`), response (`y`) and residuals (`resid`),
and the weight ``\omega_t = \sum_i w_i^2``; `synthetic(r, h; term)` gives the column of any
coefficient in the synthetic regression it uses (above). The local projection is
approximately the time-series regression of `y` on `x` weighted by ``\omega_t`` (fn 24);
the paper plots the rescaled series ``\sqrt{\omega_t}\,\hat Y_{t+h}`` and
``\sqrt{\omega_t}\,X_t`` over time (Fig. 3).

## Comparators

| the paper's label | `panellp` and `vcov` |
|---|---|
| 1W(i): clustered by unit | `WithinFE()`; `Clustered(:unit; ssc = :none)` |
| 2W: two-way clustered | `WithinFE()`; `Clustered(:unit, :time; ssc = :none, psd = false)` |
| DK98: Driscoll–Kraay, ``M = \lceil 4((T-h)/100)^{2/9} \rceil`` | `WithinFE()`; `DriscollKraaySE(:as2026; ssc = :none)` |
| t-HR: clustered by period, no lags | `WithinFE()`; `TLAHR(refine = false)` |
| t-HAR: equal-weighted cosine | `WithinFE()`; [`CosineHAR`](@ref)`()` |
| t-LAHR, plain: one lag, normal | `LagAugmented(ConstantLags(1))`; `TLAHR(refine = false)` |
| t-LAHR, recommended | `LagAugmented()`; `TLAHR()` |

The paper sets a negative two-way variance's standard error to zero (supplement D); with
`psd = false` it is NaN here, and neither covers. `CosineHAR()` uses
``B = \lfloor 0.41\,T^{2/3} \rfloor`` cosine terms and t(``B``) critical values, as
Stata's `harreg` (Lazarus, Lewis, Stock & Watson 2018); the paper says only "LLSW
tuning".

## Monte Carlo: Table 1

The stylized design of the paper's §2: ``y_{it} = 0.7 y_{i,t-1} + \beta_i X_t + Z_t +
u_{it}``, ``\beta_i \sim N(1, 1)``, ``s_i = \beta_i``, so ``\beta(h) = 0.7^h``;
``T = 30``, ``N = 1500``; ``\sigma_u^2 = N(1 - \bar R^2)/\bar R^2`` (fn 10). Coverage of
90% intervals at ``h = 0, 1, 2, 5`` over 2000 replications (the `:mc_full` test item; the
paper: 5000 replications), with the paper's values below each row:

| ``\bar R^2 = 0.9`` | ``h = 0`` | 1 | 2 | 5 |
|---|---|---|---|---|
| 1W(i) | .482 | .376 | .306 | .371 |
| *paper* | *.494* | *.372* | *.344* | *.376* |
| 2W | .812 | .726 | .726 | .804 |
| *paper* | *.819* | *.726* | *.725* | *.811* |
| DK98 | .706 | .720 | .689 | .774 |
| *paper* | *.716* | *.725* | *.695* | *.782* |
| t-LAHR plain | .832 | .830 | .801 | .780 |
| *paper* | *.851* | *.842* | *.811* | *.806* |
| t-LAHR recommended | .902 | .906 | .870 | .874 |
| *paper* | *.903* | *.899* | *.872* | *.875* |

| ``\bar R^2 = 0.1`` | ``h = 0`` | 1 | 2 | 5 |
|---|---|---|---|---|
| 1W(i) | .884 | .868 | .882 | .870 |
| *paper* | *.888* | *.878* | *.866* | *.880* |
| 2W | .824 | .812 | .808 | .804 |
| *paper* | *.841* | *.829* | *.814* | *.820* |
| DK98 | .814 | .801 | .782 | .800 |
| *paper* | *.809* | *.791* | *.789* | *.810* |
| t-LAHR plain | .829 | .828 | .828 | .802 |
| *paper* | *.848* | *.849* | *.836* | *.835* |
| t-LAHR recommended | .900 | .901 | .900 | .898 |
| *paper* | *.903* | *.902* | *.891* | *.895* |

The information criterion selects one lag in 85% (``\bar R^2 = 0.9``) and 80%
(``\bar R^2 = 0.1``) of the replications. The fn 10 calibration reproduces the paper's
unit-cluster row. By eq. (7), whose exact calibration is ``\sigma_u^2 = N(1 - \bar
R^2)(1 + \rho^2\bar\beta^2)/\bar R^2``, it corresponds to ``\bar R^2 \approx 0.93``; with
the exact calibration the unit-cluster coverage rises to .554, .440, .368, .439, further
from the paper. A reduced version (``N = 200``, 300 replications) runs in the test suite.

## Cross-checks

- **The authors' `panel_LP.R`**, run on our two synthetic panels (balanced, and with 10%
  of the rows removed): five specifications (``s_i`` with unit and time effects;
  ``s \equiv 1`` with unit effects; a time-varying exposure; a control; two exposures),
  ``p_{\max} \in \{0, 2\}``, level and cumulative responses, with and without the
  refinement: 672 rows. With `LagAugmented(CappedLags(pmax))`, `singletons = :keep` and
  `TLAHR(weights = :raw)`, they agree on the balanced panel within 2e-12 (estimates;
  standard errors 1e-13, degrees of freedom 3e-13, relative) and on the unbalanced one
  within 2e-9 (their absorption stops at a tolerance of 1e-8). The exception is the second
  exposure's degrees of freedom, where their code uses the first coefficient's contrast.
  The paper's weights agree with theirs where ``\omega_t`` is constant (above).
- **Stata.** `reghdfe` reproduces the lag-augmented two-way regression and its
  time-clustered variance (CR1, with the observation, cluster and degrees-of-freedom
  counts), within 2e-14 (the unbalanced panel's standard errors within 5e-12). On the
  synthetic regression (aggregate regressors, a balanced panel: ``s \equiv 1`` with unit
  effects, and ``\hat s_i`` with unit and time effects), `regress …, noconstant vce(hc2,
  dfadjust)` reproduces `TLAHR()`'s standard errors and degrees of freedom, and
  `harreg …, estimator(ewc)` reproduces `CosineHAR()`, within 3e-14.

| reference | `panellp` and `vcov` |
|---|---|
| `panel_LP(y, X, s, i, t, FE = cbind(id, t), p_max = p, small_sample = FALSE)` | `leads(y) ~ s & X + fe(id) + fe(t)`, `LagAugmented(CappedLags(p))`, `singletons = :keep`; `TLAHR(refine = false)` |
| `… small_sample = TRUE` | `TLAHR(weights = :raw)` |
| `… cumulative = TRUE` | `cumul(y) ~ …` |
| `reghdfe F#.y sX L(1/p).y L(1/p).sX, absorb(id t) vce(cluster t)` | `LagAugmented(CappedLags(p))`; `Clustered(:time)` |
| `regress ȳ x̄ L.x̄ L2.x̄, noconstant vce(hc2, dfadjust)` (synthetic, demeaned) | `leads(y) ~ X + lags(X, 1:2) + fe(id)`; `TLAHR()` |
| `harreg ȳ x̄ L.x̄ L2.x̄, noconstant` | the same; `CosineHAR()` |
