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
using RecipesBase: RecipesBase
using Statistics: Statistics
using StatsAPI: StatsAPI
using StatsFuns: StatsFuns
using StatsModels: StatsModels
using Tables: Tables

end # module
