# Audit: Ugarte-Ruiz (2026)

*Panel local projections without fixed-effects: a cumulative-difference OLS estimator*,
BBVA Research working paper, 29 June 2026. The paper ships no code. Its Monte Carlo runs
with the Stata commands `xtarsim` (Bruno 2005) and `locproj` (Ugarte-Ruiz 2023), its
split-panel jackknife with `locproj …, m(xtspj)` (Sun & Dhaene 2019; L369), and the
applications' jackknife with `xtlp` (Mei, Sheng & Shi; L764). Line numbers refer to the
paper converted to markdown. The findings below come from reading the paper, the
installed sources of `xtarsim` and `locproj`, and simulations in this package's test
suite.

## The estimator is as stated; the stated algebra is not

| item | paper | what holds | `panellp` |
|---|---|---|---|
| ``\gamma_h``, the coefficient on ``s_{i,t-1}`` | ``(\rho_y^{h} - 1)\beta`` (L261); ``(\rho_y^{h+1} - 1)\beta`` in eqs. 29, 34, 35 | L261 is a typo, and both are the structural moving-average coefficient, not the probability limit, which is ``(\rho_y^{h+1} - 1)\beta / (1 - \rho_y \rho_s)`` (tested) | ``\gamma_h`` is a nuisance coefficient; not reported as structural |
| Proposition 1 | "unbiased" (L269) | the proof uses population moments: it shows consistency, not finite-sample unbiasedness | documented as consistent |
| Assumption 1 | ``E[s_{i,t} \mid v_{i,t}] = 0`` (L110), ``E[s_{i,t} \mid v_{i,s}] = 0`` for ``s \ne t`` (L112) | contradicts eq. (2), since ``s_{i,t}`` contains ``v_{i,t}`` and its past; what is needed is that ``e`` is independent of the path of ``s`` | — |
| L233 | ``s_t`` uncorrelated with ``s_{t-k}`` for ``k > 2`` | false for an AR(1); the argument needs the Markov property (no partial correlation given ``s_{t-1}``) | — |

## The fixed-effects and jackknife comparators have another estimand

In the specification without the shock's lag (eq. 24; Table 9), the fixed-effects and
split-panel jackknife LPs in levels do not target ``\beta_h``: as ``T`` grows they
converge to ``\beta_h + \beta\rho_y^{h+1}\rho_s/(1 - \rho_y\rho_s)``, +177.8% of ``\beta``
at ``h = 0`` when ``\rho_y = \rho_s = 0.8`` (a test item checks this limit on a long
panel). Table 9's jackknife "bias" of +177.45% at ``h = 0`` is this omitted-lag term, not
an over-correction by the jackknife. The paper ties the upward bias to specifications
without shock lags (L834), and with the shock's lag (Table 10) the jackknife's bias at
``h = 0`` is −5.57%. The term does not explain everything: the +1 y-lag specification
(Table 11) still reports an average jackknife bias of +12.6% (+20.4% at ``h = 9``).

The Monte Carlo's jackknife (`locproj …, m(xtspj)`) is not `xtlp`'s, which this package
reproduces. `xtspj` averages two splits when a unit has an odd number of rows (here at
every odd ``h``, with ``T = 50``), and reports the OLS variance of the within regression
with normal critical values. So the paper's jackknife coverage comes from another
variance.

## Unstated choices the reported numbers depend on

- **Variance.** The paper calls the standard errors "robust" (L335–342) and recommends a
  robust or Driscoll–Kraay variance (L281), without a formula. The CD-OLS coverage
  pattern in its tables (about 0.95 with ``\rho_y = 0``, up to 1.00 at long horizons with
  ``\rho_y = 0.8``) is what heteroskedasticity-robust (HC1) standard errors give when the
  long-difference error is serially correlated. Simulating the paper's designs (a test
  item; 2000 panels for the figures): with ``\rho_s = 0``, ``\rho_y = 0.8``, HC1
  overstates the sampling standard deviation by 19–51% at ``h = 2``–``9`` (coverage
  98–99.8%); with ``\rho_s = \rho_y = 0.8`` it understates it by 7–12% at
  ``h = 1``–``4`` (coverage about 93%). Unit-clustered standard errors match it within 4%.
- **The DGP.** The signal-to-noise ratio (L305), ``\sigma_e`` and the start-up are not
  defined in the paper. They are `xtarsim`'s: ``\sigma_e = 1``, a signal-to-noise ratio
  over ``\sigma_e^2`` that sets ``\sigma_\xi`` (2.68, 2.56 and 0.72 for ``(\rho_s,
  \rho_y) = (.8, 0), (0, .8), (.8, .8)``),
  and a stationary start-up without burn-in. The test suite checks the committed
  `xtarsim` panels against these formulas, and that the package's own generator, at the
  same ``\sigma_\xi``, reproduces their moments.
- **Time effects.** Eq. (22) has a time effect ``u_t`` (L297), but the theory excludes
  them (L118) and the Monte Carlo specifications have none (L344).

## Labels and reproducibility

- Table 11 is labelled "Lags ``S_t`` 0, Lags ``y_t`` 0"; by its position it is the
  "+1 y lag" specification.
- Table A12 is labelled ``\rho_y = 0``; by its position it is the ``\rho_s = \rho_y = 0.8``
  design.
- Appendix A.5 reads "0.8 instead of 0.8" (L1164); it means 0.5.
- CMRE with one added shock lag is the same regression as CMRE1 without added lags, yet
  the tables report different numbers, so the specifications do not share simulated data
  (or the implementations differ). Exact reproduction of the tables is not possible.

## Limits of the method the paper does not test

- **Heterogeneous shock means.** The design's shock has the same mean for every unit.
  With unit-specific means the pooled regression of ``s_{i,t}`` on ``s_{i,t-1}`` keeps the
  mean in the residual, and CD-OLS is inconsistent. The test suite shows a bias of tens of
  standard errors on a large simulated panel.
- **Unit roots and trends.** With ``\rho_y = 1`` or unit-specific trend growth, the long
  difference keeps a unit term.
- **Feedback** from ``y`` to ``s`` is excluded by assumption.

## `xtarsim`

Read from the installed `xtarsim.ado` (v1.0.0):

- negative parameters are mis-handled: Stata evaluates `^` before unary minus, and the
  parameters are pasted into expressions as text, so ``-0.5^2`` is ``-0.25``; with
  ``\beta < 0`` the command stops with error 411. The paper uses non-negative parameters
  only;
- with ``\gamma\rho \ne 0`` the start-up of ``x`` and of its contribution to ``y`` use the
  same first innovation, so the first periods carry a transient term decaying like
  ``\rho^{t-2}``;
- a two-way model's time effect enters the dynamic equation as
  ``(\tau_t - \gamma\tau_{t-1})/(1-\gamma)``, not ``\tau_t``;
- `seed()` seeds the kiss32 generator (the command runs under `version 8.2`), so the draws
  are reproducible across Stata versions but not from Julia; the panels are committed.
