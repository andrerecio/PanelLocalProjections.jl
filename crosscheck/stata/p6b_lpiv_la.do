* Phase 6b cross-check: lag-augmented panel LP-IV (Almuzara & Sancibrián 2026, §3.4) on the
* fx_lpiv panels. At horizon h, p = min(h, 2) lags of y (exogenous) and of s_i ffr_t
* (endogenous, instrumented by the same lags of s_i mp_t) are added:
*     ivreghdfe Fy w L(1/p).y (sffr L(1/p).sffr = smp L(1/p).smp), absorb(…) vce
* A. unit and time effects: clustered by period, robust, Driscoll–Kraay.
* B. unit effects only, clustered by period (no fixed effect nested in the clusters).
* Writes test/fixtures/stata/p6b_lpiv_la.csv; spec labels map to panellp/vcov calls in
* test/test_crosscheck_p6.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p6b_lpiv_la.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df,N,G" _n

capture program drop _row
program define _row
    args fh ds spec h term
    local G = cond(missing(e(N_clust)), ., e(N_clust))
    local df = cond(missing(e(df_r)), ., e(df_r))
    file write `fh' "`ds',`spec',`h',`term'," %25.17e (_b[`term']) "," ///
        %25.17e (_se[`term']) "," %12.0f (`df') "," %12.0f (e(N)) "," %12.0f (`G') _n
end

local H 4
foreach ds in fx_lpiv fx_lpiv_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    generate double sffr = s * ffr
    generate double smp = s * mp
    forvalues h = 0/`H' {
        quietly generate double Fy`h' = F`h'.y
        local p = min(`h', 2)
        local ly ""
        local le ""
        local li ""
        if `p' > 0 {
            local ly "L(1/`p').y"
            local le "L(1/`p').sffr"
            local li "L(1/`p').smp"
        }
        local m "Fy`h' w `ly' (sffr `le' = smp `li')"
        quietly ivreghdfe `m', absorb(id t) cluster(t)
        _row `fh' `ds' la_twfe_cl_t `h' sffr
        quietly ivreghdfe `m', absorb(id t) robust
        _row `fh' `ds' la_twfe_robust `h' sffr
        quietly ivreghdfe `m', absorb(id t) dkraay(3)
        _row `fh' `ds' la_twfe_dk3 `h' sffr
        quietly ivreghdfe `m', absorb(id) cluster(t)
        _row `fh' `ds' la_fe_cl_t `h' sffr
    }
}
file close `fh'
display as result "p6b_lpiv_la: done"
