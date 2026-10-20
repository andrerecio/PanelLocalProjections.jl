# Split-panel jackknife

The fixed-effects (within) estimator of a panel local projection is biased in panels with a
short time dimension (Nickell 1981). In a local projection the bias is present even
without a lagged dependent variable. At ``h ≥ 1`` the error is a moving average of future
innovations, so it is not strictly exogenous, and the within transformation brings
future values of the regressors into every observation. The bias is ``O(1/T)`` and
distorts tests even when ``N`` and ``T`` grow together (Mei, Sheng & Shi 2026,
Proposition 1). For a persistent shock (``\rho > 0``) it attenuates the response and
grows with ``h`` and with the persistence (Remark 3).

[`SplitPanelJackknife`](@ref) removes the leading term:

```math
\hat\theta^{\,spj}_h = 2\hat\theta_h - \tfrac12\big(\hat\theta_{h,a} + \hat\theta_{h,b}\big)
```

``\hat\theta_h`` is the within estimate on horizon ``h``'s sample, and
``\hat\theta_{h,a}``, ``\hat\theta_{h,b}`` are the within estimates on its two halves
(Mei, Sheng & Shi 2026, eq. 16).

```julia
r = panellp(@formula(leads(y) ~ x + lags(x, 1:4) + lags(y, 1:4) + fe(id)), data;
            id = :id, time = :t, horizon = 10, estimator = SplitPanelJackknife())
summarize(r, Clustered(:unit; ssc = :none))   # the paper's unit-clustered variance
```

## The halves

The split is made **per unit, on each horizon's own estimation sample**. A unit's rows are
taken in time order:

- `split = :ceil` (default) puts the first ``\lceil n_i/2\rceil`` rows in half ``a``, as
  the authors' code and Stata `xtlp` do;
- `split = :floor` puts the first ``\lfloor n_i/2\rfloor`` there, as in the paper's
  ``\mathcal T_a = \{t \le T_h/2\}``.

The two rules differ only when ``n_i`` is odd, and the choice moves published numbers: on
the Romer–Romer GDP sample the paper's own split puts the Table 1 peak at −6.606, not
−6.285 (see the [replication page](../replication/mss2026.md)). Appendix B.2.4 splits
each unit's contiguous span of an even length. The per-unit count split, which is the
rule of the authors' code and of `xtlp`, extends this to odd counts and to gaps inside a
unit. **Every fixed effect is re-estimated inside each half**, including time effects,
which are then specific to the half (appendix B.2.3). The halves never drop rows. Under
exact absorption a unit with a single row in a half contributes nothing to that half's
estimate. Under `OneSweep()` with time effects, its demeaned zero still enters the period
means of its half, as in the authors' code. A regressor that
the fixed effects absorb inside a half (but not in the full sample) leaves that half
rank-deficient; the horizon is then `:spj_half_singular`, with `NaN` coefficients.

With fixed effects other than the unit, the singleton policy changes the halves. A row
that is a singleton in another dimension, a period observed for one unit only say, is
absorbed in the full sample but still counts in its unit's ``n_i``. So
`singletons = :drop` (as `xtlp`) and `:keep` (as the authors' R code) give different SPJ
estimates. With unit effects alone they agree.

Each horizon absorbs three samples: the full one and its two halves. With two-way effects
the thinner halves need more sweeps, so SPJ takes two to five times as long as
`WithinFE`; on the 1.1-million-row benchmark (firm and sector × quarter effects, 21
horizons) it takes about 5×.

## Inference

The variance is a sandwich with the **full-sample bread** ``(\tilde X'\tilde X)^{-1}``,
the **SPJ residuals** ``\tilde y - \tilde X\hat\theta^{spj}`` on the full-sample within
data, and the **scores** ``d = 2\tilde X - \tilde X_{\text{half}}``. Here
``\tilde X_{\text{half}}`` holds each row within-transformed inside its own half. On a
balanced panel with an even number of periods and unit effects, ``d`` is the paper's
``d^*`` (Assumption 2): ``W_{it} - \bar W_{i,b}`` in half ``a`` and
``W_{it} - \bar W_{i,a}`` in half ``b``. With time effects as well, it is the two-way
``d^*`` of appendix B.2.3, that expression minus its cross-sectional mean at ``t``. The
test suite checks both, with variances that tell the scores apart row by row. For an odd
``n_i``, ``d = W_{it} - 2\bar W_i + \bar W_{i,\text{own half}}``, as in the authors' code
and `xtlp`.

Every variance specification applies unchanged. The paper uses unit clustering without a
small-sample factor (`Clustered(:unit; ssc = :none)`). The SPJ and FE estimators have the
same asymptotic variance, but SPJ intervals can be wider in finite samples.

The authors' rule of thumb (§2.3) is to use SPJ when ``T ≥ 30`` and ``N/T ≤ 10``. Their
theory requires ``N/T^3 → 0`` and a fixed ``h``.

## Cross-checks

| reference | what is compared | agreement |
|---|---|---|
| authors' `LP_panel(method = "SPJ")` on the synthetic panels (h = 0–3): unit and two-way FE, one-way and two-way clusters, their `robust`/`eigens` options | coefficient and SE | 1e-14 |
| authors' replication notebook on their four applications (local data) | 128 coefficients and SEs | about 1e-12 (relative) |
| Stata `xtlp …, fe m(spj)` (h = 0–3, both panels): `vce(cluster id [, ase])`, `cluster t`, `cluster reg`, `cluster id t [, ase]`, `robust`, `dkraay` (default lags; `lag(2), ase`; `lag(2), nodfadj`), `keepsingletons`; without lags (single-row halves), drop and keep | coefficients, SEs, N, df, clusters | b 1e-13, se 9e-15 |
| Stata `xtlp …, tfe m(spj)`: `vce(cluster id [, ase])`, `cluster t`, `cluster id t, ase`, `robust`, `dkraay lag(2), ase`, `keepsingletons`; without lags, drop and keep | coefficients, SEs, N, df, clusters | b 1e-13, se 2e-11 (xtlp's absolute stopping rule) |
| the formula ``2\hat\theta - \tfrac12(\hat\theta_a + \hat\theta_b)`` with hand-made halves, and the variance from explicit dummies (fx_unbalanced; `:ceil` and `:floor`; unit and two-way FE; with and without lags) | coefficients; unit-clustered and robust variances | 1e-10 |
| Stata `xtspj, model(regress) method(parm)` (Sun & Dhaene 2019): fx_balanced, h = 1, 3, where every unit has the same even ``n_h`` | point estimates | 6e-14 |
| Monte Carlo of the paper's toy design (``\rho = 0.8``, ``N = 30``, ``T = 60``), 300 replications, h = 1, 5, 10 | mean IRF and 95% coverage of FE and SPJ against the authors' 1000 replications; SPJ standard error against the sampling sd | within Monte Carlo error |

`xtlp`'s `vce(unadjusted)` has no counterpart: `panellp` has no homoskedastic variance.
