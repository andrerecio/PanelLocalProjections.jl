# Audit: Ramey's infrastructure replication code

Replication files of Ramey, "The Macroeconomic Consequences of Infrastructure
Investment" (README dated June 5, 2021): Dynare models with Stata programs that plot
their impulse responses and compute multipliers (`Models/`), and Stata programs for the
empirical figures (`Empirical/`). They are a local reference here (not committed). The
findings come from reading the code; line numbers refer to the files named.

## The local projection is a cross-section LP-IV

The only local projection, `Empirical/arra_highway_irfs.do` (Figs. 6–7), is not a
pooled panel local projection. Each regression is a cross-section of states in 2009
(Alaska and Nebraska dropped, l.110). The panel only serves to build long differences:

```math
y_{i,2009+h} - y_{i,2008} = a_h + \beta_h\,\mathrm{ARRA}_{i,2009} + \gamma_h'\Delta c_{i,2009}
  + \delta_h \Delta \mathrm{HE}_{i,2008} + u_{ih}, \qquad h = 0, \dots, 4,
```

where ARRA is the highway apportionment per capita and HE highway employment. The
apportionment is instrumented by two road factors (`$road_factors`, l.32), and the
controls are political changes (l.132) and lagged highway employment (l.159). It is
estimated by two-step efficient GMM (`ivregress gmm …, robust`, l.168–175) with 90% bands
of ±1.65 standard errors (l.193–194). The outcomes are ARRA obligations and outlays,
highway spending and highway employment.

`panellp` has no instrumental-variable estimator yet (planned for v0.2), so this
specification is not reproduced.

## Findings

- **The instruments.** A comment says the standard errors use "the 3 instruments" (l.134),
  but the regressions use the two of `$road_factors`. The three-factor list
  `L1_all_factors` (l.36) is defined and never used.
- **Missing values become zeros.** `replace FHWA_<name> = 0 if FHWA_<name> == .` (l.41)
  applies to every FHWA flow. A missing later year then enters a long difference as a
  spurious 0.
- **Flows, not cumulative multipliers.** The responses are long differences of annual
  flows. The cumulative variables (`cumadd_0910_*`, `cumadd_0912_*`, l.42–45, 68) are
  built and never used.
- **A control left out.** The file notes that the control `S12ces_sa_pc1208_n` "is
  significant and affects results" (l.125), and the baseline excludes it. The
  population control of Fig. 7 is a manual edit of l.159/162.
- **Weak instruments.** Two-step GMM on 48 states, with normal critical values. The
  first-stage regression is printed (l.151), but no weak-instrument or
  overidentification test accompanies the estimates.
- **Unused results.** The population-weighted regressions (l.174–178) are computed and
  not plotted (the figure combines the unweighted ones, l.207).
- **Saved settings.** The do-file has no `xtset`. Its leads and lags rely on the panel
  settings stored in the `.dta` file.
- **`Empirical/bouakez_reg.do`** (§3.3) runs time-series regressions of TFP on capital
  stocks with leads of their differences and a trend, in the style of dynamic OLS, with
  heteroskedasticity-robust (HC1) standard errors (`vce(robust)`, l.30–34). Inference on
  a cointegrating regression needs long-run-variance standard errors. The time index is
  built from the row order (`ydate = y(1960) + _n − 1`, l.18).
- **`Models/graphirfs_nk.do`, l.55.** `replace grossrproduct = grossr*grossrproduct[_n-1]`
  refers to `grossr`, which no variable is named. Stata abbreviates it to
  `grossrproduct`, so the running product is squared and the present-value multiplier of
  that file is wrong. The figures plot impulse responses only. The Figure 4 code that the
  README assigns to this file is commented out (l.166–220).
- **The README** names `Simulations_irf.xlsx` (the file is `Simulation_irfs.xlsx`),
  `neo_iadj_util.mod` (`neo_iadj_utiliz.mod`), and `Multipliers.do`, `Graphs_bea_govk.do`
  and `Bouakez_reg.do` (lower case on disk).
- **Timing in the model files.** The ``h = 0`` row of the Dynare output is the steady
  state, and the shock hits at ``h = 1`` (sheet `neo_gc` of `Simulation_irfs.xlsx`).
  So a cumulative multiplier "at ``h = 4``" covers the impact quarter and three more.
