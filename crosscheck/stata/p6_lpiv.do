* Phase 6a cross-check: panel LP-IV (two-stage least squares on the fixed-effects design)
* on the fx_lpiv panels.
*
* A. ivreghdfe: y_i,t+h on (s_i ffr_t = s_i mp_t), w, L.y, absorbing unit and time effects;
*    robust, clustered by unit, by period, two-way, and Driscoll–Kraay; each without and
*    with `small`. Both panels.
* B. ivreghdfe absorbing unit effects only, and xtivreg2, fe (the same model; it clusters
*    only on groups nested in the panels).
* C. Several endogenous regressors: g_k ffr_t instrumented by g_k mp_t (k = 1..3), unit
*    effects, clustered by period.
* Also records the Kleibergen–Paap rk Wald F (e(widstat)) for the first-stage diagnostics.
* Writes test/fixtures/stata/p6_lpiv.csv; spec labels map to panellp/vcov calls in
* test/test_crosscheck_p6.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p6_lpiv.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df,N,G,F" _n

capture program drop _row
program define _row
    args fh ds spec h term
    local G = cond(missing(e(N_clust)), ., e(N_clust))
    local df = cond(missing(e(df_r)), ., e(df_r))
    file write `fh' "`ds',`spec',`h',`term'," %25.17e (_b[`term']) "," ///
        %25.17e (_se[`term']) "," %12.0f (`df') "," %12.0f (e(N)) "," ///
        %12.0f (`G') "," %25.17e (e(widstat)) _n
end

local H 4
foreach ds in fx_lpiv fx_lpiv_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    generate double sffr = s * ffr
    generate double smp = s * mp
    generate double Ly = L.y
    forvalues k = 1/3 {
        generate double g`k'ffr = g`k' * ffr
        generate double g`k'mp = g`k' * mp
    }
    forvalues h = 0/`H' {
        quietly generate double Fy`h' = F`h'.y
        local m "Fy`h' w Ly (sffr = smp)"
        * --- A. two-way fixed effects
        foreach sm in "" "small" {
            local tag = cond("`sm'" == "", "", "_small")
            quietly ivreghdfe `m', absorb(id t) robust `sm'
            _row `fh' `ds' twfe_robust`tag' `h' sffr
            quietly ivreghdfe `m', absorb(id t) cluster(id) `sm'
            _row `fh' `ds' twfe_cl_id`tag' `h' sffr
            quietly ivreghdfe `m', absorb(id t) cluster(t) `sm'
            _row `fh' `ds' twfe_cl_t`tag' `h' sffr
            quietly ivreghdfe `m', absorb(id t) cluster(id t) `sm'
            _row `fh' `ds' twfe_cl_id_t`tag' `h' sffr
            quietly ivreghdfe `m', absorb(id t) dkraay(3) `sm'
            _row `fh' `ds' twfe_dk3`tag' `h' sffr
        }
        * --- B. unit effects: ivreghdfe and xtivreg2, fe
        quietly ivreghdfe `m', absorb(id) cluster(t)
        _row `fh' `ds' fe_cl_t `h' sffr
        * (xtivreg2, fe clusters only on groups nested in the panels)
        quietly xtivreg2 `m', fe cluster(id)
        _row `fh' `ds' xtivreg2_fe_cl_id `h' sffr
        quietly xtivreg2 `m', fe robust small
        _row `fh' `ds' xtivreg2_fe_robust_small `h' sffr
        * --- C. several endogenous regressors
        quietly ivreghdfe Fy`h' w (g1ffr g2ffr g3ffr = g1mp g2mp g3mp), absorb(id) cluster(t)
        forvalues k = 1/3 {
            _row `fh' `ds' grp_fe_cl_t `h' g`k'ffr
        }
    }
}
file close `fh'
display as result "p6_lpiv: done"
