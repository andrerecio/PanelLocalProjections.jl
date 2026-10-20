# Inference

Variances are computed after estimation, from score sums stored at every horizon. Many
specifications can therefore be compared without re-estimating:

```julia
r = panellp(f, data; id = :id, time = :t, horizon = 10)
summarize(r, Clustered(:unit))                       # one-way, by unit
summarize(r, Clustered(:unit, :time))                # two-way (Cameron–Gelbach–Miller)
summarize(r, DriscollKraaySE())                      # Driscoll–Kraay
cv = vcov(Robust(), r); vcov(cv, 3)                  # the matrix at horizon 3
```

## Specifications

| specification | variance |
|---|---|
| `Robust(; ssc)` | heteroskedasticity-robust |
| `Clustered(dims...; ssc, psd)` | one- or two-way cluster-robust. Dimensions are `:unit`, `:time`, a fixed-effect or declared cluster column, or a tuple of columns for their intersection. Columns other than the unit, the time and the fixed effects must be declared at estimation with `clusters = (:region,)` |
| `DriscollKraaySE(bw; ssc, timelags)` | Driscoll–Kraay: a Bartlett long-run variance of the per-period score sums, ordered by calendar period |

## Small-sample conventions

The reference implementations disagree on the small-sample factor and the critical values.
Each convention is an explicit `ssc` value. Below, ``n`` is the number of observations,
``K`` the number of estimated coefficients, ``df_a`` the absorbed degrees of freedom,
``G`` the (smallest) number of clusters and ``T`` the number of periods in the sample.

| spec, `ssc` | factor | critical values | reproduces |
|---|---|---|---|
| `Robust`, `:none` | 1 | normal | HC0, fixest `ssc(K.adj = FALSE)` |
| `Robust`, `:stata` | ``n/(n-K-df_a)`` | ``t(n-K-df_a)`` | `regress`/`reghdfe`/`xtlp`, `vce(robust)` |
| `Clustered`, `:none` | 1 | normal | the papers' unit-clustered variances; the variance (not the critical values) of `xtlp …, ase`; fixest `ssc(K.adj = FALSE, G.adj = FALSE)` |
| `Clustered`, `:stata` | ``\frac{n-1}{n-a-K-df_a^{c}}\frac{G}{G-1}`` | ``t(\min(n-K-df_a^{c},\,G-1))`` | `reghdfe`, `xtlp`, `xtreg, fe vce(cluster)`, `regress, vce(cluster)`, fixest's default (one-way) |
| `Clustered(:unit)`, `:mss` | ``\frac{N}{N-1}\frac{n-1}{n-K}`` | normal | Mei–Sheng–Shi code, `robust = T` |
| `Clustered(:unit, :time)`, `:mss` | ``\frac{G}{G-1}\frac{n-1}{n-N-K}``, ``G = \min(N, T)`` | normal | Mei–Sheng–Shi code, `robust = T` |
| `DriscollKraaySE`, `:none` | 1 | normal | Mei–Sheng–Shi appendix B.2.2; the variance (not the critical values) of `xtlp`/`xtscc …, ase` |
| `DriscollKraaySE`, `:stata` | ``\frac{n-1}{n-K-df_a^{dk}}\frac{T}{T-1}`` | ``t(\min(n-K-df_a^{dk},\,T-1))`` | `xtlp`, `reghdfe vce(dkraay)` |
| `DriscollKraaySE`, `:xtscc` | ``\frac{n-1}{n-K-1}\frac{T}{T-1}`` (with fixed effects) | ``t(T-1)`` | `xtscc`, `xtlp …, nodfadj` |

The two `:mss` rows use ``N`` and ``T``, the numbers of units and periods in the data;
`ssc = :mss` exists only for these two specifications, as in the authors' code. With
unbalanced data the two-way factor turns negative at long horizons (``n \le N + K``),
which is reported. `ase` in `xtlp` and `xtscc` removes the small-sample factor but keeps
their ``t(G-1)`` or ``t(T-1)`` critical values. `summarize(…; dist = :normal)` gives
normal ones.

A variance is undefined at a horizon with a single cluster or period, or without
residual degrees of freedom. It is then `NaN` at that horizon, and the horizons are
reported.

Fixed effects nested in a cluster dimension (for example unit effects with unit clusters)
are not counted in ``df_a^{c}``, as in reghdfe. If any is nested, ``a = 1`` and the first
remaining fixed effect is charged one level. ``df_a^{dk}`` treats time as the cluster
variable, without the ``a`` term, as `reghdfe vce(dkraay)` and `xtlp` do. With three or
more fixed effects, ``df_a`` follows reghdfe's pairwise rule, which is a conservative upper
bound.

## Driscoll–Kraay bandwidth

`bw` sets the number of lags ``m`` (weights ``1 - j/(m+1)``):

| `bw` | lags | source |
|---|---|---|
| `:nw94` | ``\lfloor 4(T/100)^{2/9} \rfloor`` | `xtscc`, `xtlp` |
| `:as2026` | ``\lceil 4(T/100)^{2/9} \rceil - 1`` | Almuzara & Sancibrián's kernel bandwidth ``M`` (weights ``1 - \ell/M``) |
| `:mss2026` | ``\lfloor T^{1/4} \rfloor`` | Mei, Sheng & Shi, appendix B.2.2 |
| integer | as given | `reghdfe vce(dkraay bw)` uses `bw` = lags + 1 |

`timelags = :calendar` (the default) lags by calendar period, so a period absent from the
sample contributes a zero sum. `timelags = :position` lags over the periods present, as
`xtlp`, `xtscc` and `reghdfe` do.

## Two-way variances that are not positive semi-definite

`psd` repairs a two-way variance by clipping negative eigenvalues to zero (Cameron,
Gelbach & Miller 2011):

| `psd` | clipping | as in |
|---|---|---|
| `:stata` (default for two dimensions) | on ``S V S``, with ``S`` the diagonal of the raw regressors' standard deviations, floored at 1e-3. With an intercept, the slope block is clipped on its own and the intercept enters with scale 1 | `reghdfe` and `xtlp` |
| `:raw` | on ``V`` itself | the Mei–Sheng–Shi code, `eigens = T`; fixest's default `vcov_fix = TRUE` (which floors eigenvalues at 1e-16) |
| `:none` | none | fixest with `vcov_fix = FALSE` |

Without a repair, a two-way variance can have negative diagonal elements. Their standard
errors are `NaN`, and the horizons are reported. To reproduce fixest's default two-way
clustering, use `psd = :raw`.

## Critical values

`summarize` uses each specification's critical values: ``t`` with the degrees of freedom
above, or normal. `dist = :normal` forces normal ones, e.g. to reproduce the ±1.96 bands
of `xtlp`'s `e(irf)` or of the Mei–Sheng–Shi figures.
