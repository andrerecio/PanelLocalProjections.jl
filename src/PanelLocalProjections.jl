"""
    PanelLocalProjections

Panel local projections (Jordà 2005): for each horizon `h`, a panel regression of
`y_{i,t+h}` on a shock, with estimators and inference from

- Mei, Sheng & Shi (2026, *Journal of International Economics*), "Nickell bias in panel
  local projection": the split-panel jackknife (SPJ) correction of the fixed-effects
  estimator;
- Ugarte-Ruiz (2026, BBVA Research), "Panel local projections without fixed-effects":
  the cumulative-difference OLS (CD-OLS) estimator;
- Almuzara & Sancibrián (2026), "Micro responses to macro shocks": t-LAHR inference
  (time-clustered, lag-augmented, heteroskedasticity-robust; HC2 + Imbens–Kolesár df).

The fixed-effect absorption, the OLS solver and every variance estimator are implemented
in this package. Their design takes inspiration from Regress.jl and CovarianceMatrices.jl
(Giuseppe Ragusa), FixedEffectModels.jl (Matthieu Gomez), LocalProjections.jl (gragusa;
andrerecio fork), and the Stata commands `reghdfe` (Sergio Correia) and `xtlp` (Mei, Sheng
& Shi). None of them is a dependency.
"""
module PanelLocalProjections

using Dates: Dates
using LinearAlgebra: LinearAlgebra
using Printf: Printf
using RecipesBase: RecipesBase, @recipe, @series
using Statistics: Statistics
using StatsAPI: StatsAPI, coef, coefnames, vcov, stderror, nobs
using StatsFuns: StatsFuns
using StatsModels: StatsModels, @formula
using Tables: Tables

export panellp, PanelLocalProjection, PanelLPVcov, IRFSummary, @formula
export leads, ldiff, anchor, cumul, lags, firstdiff, fe
export WithinFE, SplitPanelJackknife, CDOLS, AlternatingProjections, OneSweep
export Robust, Clustered, DriscollKraaySE
export summarize, coefpath, nobspath, horizons, esample
export coef, coefnames, vcov, stderror, nobs
export irfplot, irfplot!

include("utils.jl")
include("types.jl")
include("panelindex.jl")
include("formula.jl")
include("design.jl")
include("absorb.jl")
include("ols.jl")
include("scores.jl")
include("estimators.jl")
include("panellp.jl")
include("vcov.jl")
include("summary.jl")
include("plotting.jl")

end # module
