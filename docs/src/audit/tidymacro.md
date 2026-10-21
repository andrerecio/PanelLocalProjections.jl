# Audit: tidyMacro `fLPPanel`

tidyMacro (Muhsin Ciftci; MIT licence; version 1.1.0, commit `dbf16c9`, 19 October 2026)
is an R package whose `fLPPanel` estimates the lag-augmented panel local projection of
Almuzara & Sancibrián with clustered or Imbens–Kolesár inference. Its C++ engine
describes its fixed-effect step as matching "`regress_HDFE` in the MATLAB / R source"
(`src/fLPPanel.cpp`, l.35), the code of Almuzara & Sancibrián, which has no licence.
This package neither uses nor runs tidyMacro. The findings below come from reading its
source at that commit; line numbers refer to `src/fLPPanel.cpp` and `R/fLPPanel.R`.

## What it computes

At each horizon, `fLPPanel` regresses ``y_{i,t+h}`` (or ``\sum_{r=0}^{h} y_{i,t+r}``
with `cumulative = TRUE`) on the interactions of the exposures with the aggregate
shocks, ``\min(h, p_{\max})`` lags of ``y`` and of the interactions, and the controls,
after absorbing the fixed effects. This is the authors' `panel_LP`, with the lag rule
`CappedLags(p_max)` (no lag at ``h = 0``, no information criterion). The variance is
clustered by period, by unit or both, with `ssc = "none"`, `"CR1"` or `"IK"` (the HC2 and
Imbens–Kolesár refinement, by period only).

## Findings

- **The degrees of freedom use the first coefficient's contrast for every coefficient**
  (`XX0_col0 = XtXt_pinv.col(0)`, l.370), as the authors' code does (see the
  [audit of that code](as_code.md), bug B2). With several shocks, every degrees of
  freedom but the first is computed for the wrong contrast.
- **The synthetic weights are the raw exposures** (l.336–347), as in the authors' code:
  ``\bar X_t = \sum_{i \in t} s_i \ddot X_{it} / \sum_{i \in t} s_i^2``, with the periods
  unweighted. The paper uses the demeaned exposure, and `TLAHR()` weights the periods by
  ``\omega_t`` (fn 24). The two differ on unbalanced panels, with time-varying exposures
  and with several exposures (see [t-LAHR](../methods/tlahr.md)). The exposure column of
  coefficient `is` is `is % n_s_cols` (l.336), which fixes the out-of-range index of the
  authors' code with several aggregate shocks (bug B1).
- **The leverage comes from a pseudo-inverse of the Gram matrix** ``\bar X'\bar X``
  (l.357–359), with Armadillo's default tolerance. When ``\bar X`` is ill-conditioned this
  squares its condition number, and the dropped directions change the leverage and the
  degrees of freedom. The residual-maker diagonal enters as ``\sqrt{|d_t|}``, floored at
  1e-14 (l.362–364), so a leverage of one does not stop the computation.
- **`Date` time columns lose most periods** (`R/fLPPanel.R`, l.47–61). Dates are converted
  to day counts and the grid step is the smallest gap between distinct dates. Rows more
  than 0.01 of a step off the grid are dropped. With quarterly dates the step is 90 days,
  so a quarter of 91 or 92 days sits 0.011–0.022 steps off; monthly dates (a 28-day step)
  fare worse. Integer periods are not affected.
- **The fixed-effect demeaning** stops when the norm of a sweep's change is below 1e-8
  times the Frobenius norm of ``[y\ X]`` (l.108–120), after at most 10,000 sweeps, without
  acceleration. Non-convergence is reported in a status code.
- **Singletons are kept**, as in the authors' code.
