# Panel LP-IV

When the macro shock of interest is not observed, but an endogenous aggregate ``\tilde X_t``
(a policy rate) and a noisy proxy ``X_t^* = a_0 X_t + \nu_t`` (high-frequency surprises)
are, Almuzara & Sancibrián (2026, §3.4) instrument the former with the latter:

```math
y_{i,t+h} = \beta^{IV}(h)\, s_i \tilde X_t + \psi(h)' W_{it} + \mu_i(h) + \nu_t(h)
            + \xi_{it}(h),\qquad s_i \tilde X_t \text{ instrumented by } s_i X_t^*,
```

with the normalisation ``b_0 = 1`` (the impact response of ``\tilde X_t`` to ``X_t``). Under
heterogeneous responses the estimand is the same as in the observed-shock case: the first
stage is aggregate, so the instrument corrects measurement error without reweighting the
units (their Proposition 3).

## The fit

A block `(endo ~ inst)` in the formula makes every horizon's fit two-stage least squares
on its fixed-effects design:

```julia
r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), data;
            id = :id, time = :t, horizon = 12)
```

- The endogenous regressors are on the left of the block, the excluded instruments on its
  right; the other regressors are exogenous and instrument themselves. Several endogenous
  regressors are allowed (e.g. group dummies × an aggregate). They are listed first among
  the coefficients and are the default `shock`.
- The instruments are absorbed with the rest of the design. The first stage projects the
  endogenous columns on the exogenous ones and the instruments, and the fitted values X̂
  replace them; the residual is the structural one, ``\tilde y - \tilde X \hat\beta``.
  The fit stores the bread ``(\hat X'\hat X)^{-1}`` and the scores ``\hat x_{it}\hat\xi_{it}``,
  so every variance (`Robust`, `Clustered`, `DriscollKraaySE`, `CosineHAR`, `TLAHR`)
  applies to the IV scores unchanged. In the just-identified case
  ``(\hat X'\hat X)^{-1}\hat X' = (Z'X)^{-1}Z'``, which is their supplement A.4.
- An instrument absorbed by the fixed effects (e.g. an aggregate series with `fe(t)`), or
  collinear with the exogenous regressors or an earlier instrument, is left out of the
  first stage, with a note. A horizon whose instruments do not identify every endogenous
  regressor (too few left, a first-stage partial R² of round-off size, or collinear fitted
  values) is `:underidentified`, with NaN coefficients.
- IV works with `WithinFE` and `LagAugmented`. The split-panel jackknife and CD-OLS have no
  IV version.

## t-LAHR for LP-IV

[`LagAugmented`](@ref) adds, at horizon ``h``, ``p(h)`` lags of the outcome (exogenous) and
of the endogenous shock (endogenous, instrumented by the same lags of every instrument):
their ``W_{it} = (s_i\tilde X_{t-1}, y_{i,t-1}, \dots)`` and ``W^*_{it} = (s_i X^*_{t-1},
y_{i,t-1}, \dots)``. Lag ``\ell`` of the instruments enters only with lag ``\ell`` of an
endogenous shock, so the fit at ``h`` depends only on ``p(h)``.

`TLAHR(refine = false)` is their LP-IV variance: CR0 by period on the IV scores with normal
critical values (supplement A.4). The HC2/Imbens–Kolesár refinement (supplement B.2) is
defined for OLS only, so `TLAHR(refine = true)`, the default, is an error on IV fits.

With [`ICLags`](@ref), ``\hat p`` is chosen on the reduced form: the regression of the outcome
on the instruments in place of the endogenous regressors, as an OLS panel LP with the
shock's instrument as its shock (the instrument with the endogenous shock's exposure,
`s & mp` for `s & ffr`), whose exposure weights the criterion. This is an extension: the
paper's criterion (supplement B.1) is stated for an observed shock.

## First-stage diagnostics

[`firststage`](@ref) reports per horizon, with any variance specification (default
`TLAHR(refine = false)`):

- the Kleibergen–Paap (2006) rk Wald F, from the singular value decomposition of the
  normalised first stage, for one or several endogenous regressors;
- the Montiel Olea–Pflueger (2013) effective F, for one endogenous regressor;
- the first-stage coefficients of the instruments, net of the fixed effects and the
  exogenous regressors, with their standard errors.

The small-sample factor uses the first stage's degrees of freedom (exogenous regressors
plus instruments), as ivreg2, ivreghdfe and weakivtest do; `firststage` reproduces
ivreghdfe's `widstat` and weakivtest's `F_eff` (see [Cross-checks](../crosscheck.md)).

The variance matters for the first stage as much as for the second. With an aggregate
instrument the effective sample size is ``T``: clustering by unit treats the ``T`` values of
the instrument as replicated across units and can report an F in the hundreds of
thousands, where the F clustered by period reveals a weak instrument (the
[US state example](../examples/us_states.md)).

```julia
fs = firststage(r)                          # CR0 by period
firststage(r; vcov = Clustered(:unit))      # for comparison
fs.kpF, fs.effF                             # per horizon
```
