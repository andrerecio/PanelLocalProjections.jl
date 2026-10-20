* Phase 3 fixtures: panels drawn by xtarsim (Bruno 2005), the generator of the Monte
* Carlo in Ugarte-Ruiz (2026). Writes test/fixtures/stata/xtarsim_<design>.csv.
*
*   y_it = ρ_y y_i,t−1 + β s_it + η_i + ε_it,   s_it = ρ_s s_i,t−1 + ξ_it,   β = 1
*   η_i = 10 (1 − ρ_y)(1 + s̄_i − s̄)   (oneway(corr 10): the paper's correlated effects)
*   signal-to-noise ratio 20, σ_ε = 1, N = 100, T = 50, stationary start-up, no burn-in
*
* Designs (ρ_s, ρ_y): s8y0 = (.8, 0), s0y8 = (0, .8), s8y8 = (.8, .8), the paper's main
* grid. xtarsim runs under version 8.2, so seed() seeds kiss32 and the draws are
* reproducible; Julia cannot regenerate them, which is why the panels are committed.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
foreach d in s8y0 s0y8 s8y8 {
    if "`d'" == "s8y0" local args "gamma(0) rho(.8) seed(20261101)"
    if "`d'" == "s0y8" local args "gamma(.8) rho(0) seed(20261102)"
    if "`d'" == "s8y8" local args "gamma(.8) rho(.8) seed(20261103)"
    xtarsim y s eta, nid(100) time(50) beta(1) snratio(20) oneway(corr 10) `args'
    rename (ivar tvar) (id t)
    format y s eta %21.0g
    export delimited id t y s eta using ///
        "$PLP_ROOT/test/fixtures/stata/xtarsim_`d'.csv", datafmt replace
}
display as result "p3_xtarsim: done"
