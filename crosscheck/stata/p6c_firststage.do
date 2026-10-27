* Phase 6c cross-check: first-stage diagnostics of panel LP-IV on the fx_lpiv panels,
* over-identified (two instruments, s_i mp_t and sv_it mp_t, for s_i ffr_t), so that the
* Montiel Olea–Pflueger effective F differs from the Kleibergen–Paap F.
* A. ivreghdfe … absorb(id t), robust / cluster(t): e(widstat).
* B. ivreghdfe … absorb(id), cluster(id) and robust: e(widstat), then weakivtest (which
*    supports ivreg2/ivreghdfe, not xtivreg2): r(F_eff).
* C. Two endogenous regressors (s_i ffr_t, sv_it ffr_t), three instruments (s mp, sv mp,
*    w mp): e(widstat), robust / cluster(t) / cluster(id).
* D. Three group dummies × ffr instrumented by × mp, no exogenous regressor (tied singular
*    values of the normalised first stage): e(widstat), robust.
* Writes test/fixtures/stata/p6c_firststage.csv; spec labels map to firststage calls in
* test/test_crosscheck_p6.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p6c_firststage.csv", write replace
file write `fh' "dataset,spec,h,stat,value,N,G" _n

capture program drop _row
program define _row
    args fh ds spec h stat v
    local G = cond(missing(e(N_clust)), ., e(N_clust))
    file write `fh' "`ds',`spec',`h',`stat'," %25.17e (`v') "," %12.0f (e(N)) "," %12.0f (`G') _n
end

local H 3
foreach ds in fx_lpiv fx_lpiv_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    generate double sffr = s * ffr
    generate double smp = s * mp
    generate double svmp = sv * mp
    generate double svffr = sv * ffr
    generate double wmp = w * mp
    generate double Ly = L.y
    forvalues k = 1/3 {
        generate double g`k'ffr = g`k' * ffr
        generate double g`k'mp = g`k' * mp
    }
    forvalues h = 0/`H' {
        quietly generate double Fy`h' = F`h'.y
        local m "Fy`h' w Ly (sffr = smp svmp)"
        quietly ivreghdfe `m', absorb(id t) robust
        _row `fh' `ds' twfe_robust `h' widstat e(widstat)
        quietly ivreghdfe `m', absorb(id t) cluster(t)
        _row `fh' `ds' twfe_cl_t `h' widstat e(widstat)
        foreach v in cluster(id) robust {
            local tag = cond("`v'" == "robust", "robust", "cl_id")
            quietly ivreghdfe `m', absorb(id) `v'
            _row `fh' `ds' fe_`tag' `h' widstat e(widstat)
            quietly weakivtest
            scalar __feff = r(F_eff)
            _row `fh' `ds' fe_`tag' `h' F_eff __feff
        }
        * --- C. two endogenous regressors, three instruments (m = 2, L = 3)
        local m2 "Fy`h' w Ly (sffr svffr = smp svmp wmp)"
        quietly ivreghdfe `m2', absorb(id t) robust
        _row `fh' `ds' m2_twfe_robust `h' widstat e(widstat)
        quietly ivreghdfe `m2', absorb(id t) cluster(t)
        _row `fh' `ds' m2_twfe_cl_t `h' widstat e(widstat)
        quietly ivreghdfe `m2', absorb(id) cluster(id)
        _row `fh' `ds' m2_fe_cl_id `h' widstat e(widstat)
        * --- D. group dummies × the aggregate, no exogenous regressor: Π̂ is diagonal
        *    with tied singular values (the SVD's trailing block is singular)
        quietly ivreghdfe Fy`h' (g1ffr g2ffr g3ffr = g1mp g2mp g3mp), absorb(id) robust
        _row `fh' `ds' diag_fe_robust `h' widstat e(widstat)
    }
}
file close `fh'
display as result "p6c_firststage: done"
