# Audit: `xtlp`

`xtlp` (Mei, Sheng & Shi; MIT licence) is the Stata command of the Mei–Sheng–Shi paper:
the panel local projection with fixed effects (`m(fe)`) or the split-panel jackknife
(`m(spj)`), with one-way (`fe`) or two-way (`tfe`) effects. The audited version is
1.1.1 (`xtlp.ado`, 1995 lines; GitHub commit `cd3fa24`). The cross-checks of this
package run it on the synthetic panels (`crosscheck/stata/p1_fe.do`, `p2_spj.do`); see
the [cross-check map](https://github.com/andrerecio/PanelLocalProjections.jl/blob/main/crosscheck/README.md).

Statements marked *verified* are matched by those cross-checks. The others come from
reading the installed source, and line numbers refer to `xtlp.ado`.

## What it computes

At each horizon, `xtlp` builds the transformed outcome (`ytransf`: level, `diff`,
`cmltdiff` y_{t+h} − y_{t−1}, `cmltsum`), restricts the sample to complete cases, drops
singletons iteratively, demeans the outcome and the regressors, and estimates by OLS
without an intercept. With `m(spj)` it splits each unit's rows in two halves of
``\lceil T_i/2 \rceil`` and ``\lfloor T_i/2 \rfloor`` (by count, after singleton
removal), re-demeans each half, and reports ``2\hat\beta - (\hat\beta_a +
\hat\beta_b)/2``.

*Verified:* coefficients, standard errors, observation and cluster counts and critical
degrees of freedom agree with `panellp` (`WithinFE`, `SplitPanelJackknife(split =
:ceil)`) within 2e-10 for every `fe`/`tfe` × `vce()` specification of the cross-checks,
including unbalanced panels, gaps, singleton cascades, and SPJ halves with single rows.

## Conventions, and how `panellp` reproduces them

| `xtlp` | `panellp` and `vcov` |
|---|---|
| singletons dropped iteratively (units; also periods under `tfe`) | `singletons = :drop` (the default); `keepsingletons` is `:keep` |
| `tfe` demeaning: unit step, then time step, stop at an absolute step below 1e-9 or after 1000 sweeps | `AlternatingProjections()` converges to the exact within transform; the difference is within the tolerance on the test panels |
| `vce(cluster id)`, `vce(cluster t)`, `vce(cluster id t)` | `Clustered(:unit)`, `Clustered(:time)`, `Clustered(:unit, :time)`; `ssc = :stata` reproduces the factor ``(n-1)/(n-a-K-df_a) \cdot G/(G-1)`` with the nested-effects rule |
| `vce(robust)` | `Robust()` |
| `vce(dkraay)`, `lag(m)` | `DriscollKraaySE(:nw94 or m; timelags = :position)`: lags count positions among the sample's periods, not calendar distance |
| `vce(dkraay …, nodfadj)` | `ssc = :xtscc` |
| `ase` | `ssc = :none` (`df_r` is unchanged) |
| two-way variance not positive semi-definite | `psd = :stata`: negative eigenvalues clipped in the metric of the regressors' standard deviations (floored at 1e-3) |
| `e(irf)` bands: ±1.96 standard errors for every `vce()` | `summarize(…; dist = :normal)` |

The cross-check do-files set `type double`: otherwise `xtlp` stores the transformed
outcomes as floats (below).

## Findings

**Precision and stored results**

- The transformed outcomes (`level`, `diff`, `cmltdiff`) are created with `generate`
  and stored in Stata's default type, float, so they carry about 7 significant digits
  (l.405–428). Only `cmltsum` is double. `set type double` restores full precision.
- `e(irf)` is filled from float variables too, so its values are float-rounded
  (l.673–676). `e(b#)` and `e(V#)` are full double.
- The `e(irf)` bands use ±1.96 for every `vce()`, while `e(df_r#)` holds t degrees of
  freedom (l.723–724).
- Runs over several horizons store no per-horizon observation or cluster count
  (`e(N#)`, `e(N_g#)`, `e(N_clust#)`). They also never post or clear `e()`, so results
  of an earlier command persist (l.708–715). A horizon's count needs a `hor(0)` run on a
  pre-built lead.
- `hor(0)`, the default, ignores `ytransf()` silently: the transforms are built only
  when `hor > 0` (l.397).

**Samples and clusters**

- Cluster variables are never marked out of the sample: rows with a missing cluster value
  stay in it and form their own cluster (`reghdfe` would drop them).
- The two-way cluster key is `C1 · (max(C2) + 1) + C2` (l.1636). It is unique only for
  non-negative integer `C2`. Negative codes, such as `%tq` or `%tm` dates before 1960,
  can merge distinct cells: with max(C2) = 100, (C1 = 2, C2 = −40) and (C1 = 1,
  C2 = 61) both map to 162.
- The Driscoll–Kraay gap warning assumes a unit spacing of the time variable and ignores
  the `xtset` delta (l.1786). With a delta of 5, e.g. data every 5 years, it warns about
  gaps that do not exist.

**Estimation**

- The `tfe` demeaning stops when no element changes by more than 1e-9 in absolute value
  (l.1936–1989). The rule depends on the scale of the data: with large values, rounding
  noise can keep it running for all 1000 sweeps. On slowly mixing unbalanced panels, the
  remaining error can be well above 1e-9. When the limit is reached, it returns without
  a warning.
- No collinear regressor is dropped. A full set of dummies, or a base level of a factor
  variable, makes ``X'X`` singular: error 198 at ``h = 0`` and a missing horizon
  otherwise.
- When a horizon's variance is undefined (no residual degrees of freedom, a single
  cluster or period), its coefficients are posted as missing too, although they could be
  computed. `panellp` keeps them and reports a NaN variance.

**Degrees of freedom**

- *Verified:* the Driscoll–Kraay default factor under `tfe` uses ``n - K - N_g + 1`` in
  its denominator, while clustering by period gives ``n - K - N_g``: the two differ even
  with no lags. `DriscollKraaySE(ssc = :stata)` reproduces both.
- The default variance (`vce(unadjusted)`) under `tfe` counts ``N_g + G_t - 1`` absorbed
  effects, which assumes one connected component of units and periods (`reghdfe` counts
  the components). `panellp` has no homoskedastic variance.
- *Verified:* inside the SPJ halves, singletons are kept, and their scores are
  ``2\tilde X_{it}``. A unit with two or three rows has a single row in a half. The halves
  are not required to have four periods per unit, although the error message for a
  singular half says so.

**Minor**

- Comments on the `ase` factors (l.1494–1691) describe factors that the code no longer
  applies (it sets them to 1).
- In the several-shocks path, the PSD-repair warning reports a stale horizon (`XTLP_H`).
