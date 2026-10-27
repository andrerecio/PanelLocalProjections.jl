# Audit: Mei–Sheng–Shi replication code

**Source.** [metricshilab/panel-lp-replication](https://github.com/metricshilab/panel-lp-replication)
at commit `882b046`, MIT licence:

- `applications/LP_panel_all.r` (the function behind every empirical result);
- `simulations/`, `simulations_appendix/`;
- `applications/replication.ipynb`.

Line numbers refer to that commit. The findings come from reading the code. Those marked
*verified* were checked by running it on the synthetic panels or on the authors' data.

## What the published numbers depend on

The notebook calls `LP_panel` with `H = 0` and outcome, lead and lag columns pre-built in
Stata. Every published empirical number comes from that mode. In it, the function is
reproduced exactly by `panellp` (*verified*, about 1e-12 on all four applications). Four
choices of the code move the published numbers:

| choice | effect (*verified*) | `panellp` |
|---|---|---|
| **One-sweep two-way demeaning** (`te = T`, lines 215–301) | exact only on a balanced sample; moves the Romer–Romer unemployment SPJ path by up to 6% and its standard errors by 12–18% (see the [replication page](../replication/mss2026.md)) | `OneSweep()` reproduces it; the default `AlternatingProjections()` is exact |
| **The split ``\lceil n_i/2\rceil``**, where the paper's ``\{t \le T_h/2\}`` gives ``\lfloor n_i/2\rfloor`` | with odd counts: Romer–Romer GDP peak (Table 1) −6.285 instead of −6.606 under the paper's split; Cerra–Saxena ``h = 10``: −7.893 vs −7.698; Baron–Verner–Xiong bank crash, ``h = 4``: −4.213 vs −4.297 | `split = :ceil` (default, as the code and `xtlp`) or `:floor` (the paper) |
| **The eigenvalue repair** (`eigens = T`) of the Baron–Verner–Xiong two-way variance, absent from the paper's formula (online appendix B.2.1) | fires at every horizon; standard errors up to 13% larger than unrepaired | `psd = :raw` |
| **File-level small-sample factors** (`robust = T`, lines 335 and 368): ``N`` and ``T_0`` count the units and periods in the input file, not in the estimation sample; the two-way factor subtracts ``N`` | changes the Baron–Verner–Xiong standard errors only | `ssc = :mss` |

## Latent bugs (not triggered by the published results)

| ID | location | problem |
|---|---|---|
| B1 | 160–162 | with `H > 0`, controls are dated at the outcome period ``t+h``, not ``t`` |
| B2 | 351–353, 964–966 | two-way clustering indexes past the end of the sample unless `H = 0`, `lagY = 0`, `lagX = 1` |
| B3 | 251, 475, 661, 836 | `te = T` indexes past the end with `lagX ≥ 2` or `lagY ≥ 2` |
| B4 | 379, 994 | `eigens <- eigen(...)` overwrites the logical argument, so `eigens = T` works for one horizon per call |
| B7 | 36–40 | periods are recoded by rank, so internally built leads and lags jump over periods missing for every unit |
| B8 | 49–58 | a unit is sorted by time only when it misses one of the file's periods; units observed in every period of the file are not re-sorted, so their rows must already be in time order |
| B9 | 72–75 | `c.name = NULL` silently uses every other column as a control |

## Discrepancies with the paper

- **Odd split.** The code puts ``\lceil n_i/2\rceil`` rows in the first half; the paper's
  ``\{t \le T_h/2\}`` puts ``\lfloor n_i/2\rfloor``. `SplitPanelJackknife(split = :ceil)`
  follows the code (and `xtlp`); `:floor` follows the paper. The effect on the published
  numbers is in the table above.
- **Small-sample factors.** The paper's variances have none. The main-text simulation
  code always applies ``N/(N-1)\,(n-1)/(n-K)``; the appendix code applies none.
- **Appendix Driscoll–Kraay.** The December 2025 fix was applied to the FE branch only.
  The SPJ branch (`simulations_appendix/scripts/02-03-error-calculation.R`, lines 158 and
  175) still divides the long-run variance by ``T_h`` and adds the unit-clustered meat,
  which affects Figures C4–C6. `panellp` implements appendix B.2.2 as written.
- **Analytical bias correction (DB).** The simulation code estimates ``\rho`` by pooled OLS,
  drops the third term of ``f_{T,h}``, and uses the FE standard error. It is not part of
  this package.

## Not used

The EPS figures were not used as targets. The plotted paths come from the same R code that
`panellp` reproduces to about 1e-12 on the same data.
