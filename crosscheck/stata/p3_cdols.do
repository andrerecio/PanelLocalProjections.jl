* Phase 3 cross-check: cumulative-difference OLS (Ugarte-Ruiz 2026),
*     y_i,t+h − y_i,t−1 = c + β_h s_it + γ_h s_i,t−1 (+ lags) + error,   pooled OLS,
* on the xtarsim panels (crosscheck/stata/p3_xtarsim.do), the native fx_cdols panel, and
* the unbalanced fx_unbalanced (shock x; gaps, singleton units, a period absent for all
* units, so no xtscc there).
* Writes test/fixtures/stata/p3_cdols.csv (one row per dataset × spec × horizon × term).
* Each spec label maps to a panellp/vcov call in test/test_crosscheck_p3.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p3_cdols.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df_r,N,G" _n

capture program drop _emit
program define _emit
    * Coefficients and standard errors of the last estimation (e(b), e(V)).
    args fh ds spec h
    foreach v in s L.s L2.s LD.y {
        capture scalar __b = _b[`v']
        if _rc continue
        scalar __se = _se[`v']
        file write `fh' "`ds',`spec',`h',`v'," %25.17e (__b) "," %25.17e (__se) "," ///
            %12.0f (e(df_r)) "," %12.0f (e(N)) "," %12.0f (e(N_clust)) _n
    }
end

capture program drop _emit_irf
program define _emit_irf
    * locproj: e(irf) has one row per horizon 0..H (IRF, Std.Err., bands); the number of
    * observations of each horizon is column 1 of the global matrix stats.
    args fh ds spec H
    matrix __irf = e(irf)
    forvalues h = 0/`H' {
        file write `fh' "`ds',`spec',`h',s," %25.17e (__irf[`h' + 1, 1]) "," ///
            %25.17e (__irf[`h' + 1, 2]) ",.," %12.0f (stats[`h' + 1, 1]) ",." _n
    }
end

capture program drop _emit_hor
program define _emit_hor
    * xtlp hor(0 H): e(b#), e(V#), e(df_r#) per horizon.
    args fh ds spec H
    forvalues h = 0/`H' {
        matrix __b = e(b`h')
        matrix __V = e(V`h')
        local j = colnumb(__b, "s")
        file write `fh' "`ds',`spec',`h',s," %25.17e (__b[1, `j']) "," ///
            %25.17e (sqrt(__V[`j', `j'])) "," %12.0f (e(df_r`h')) ",.,." _n
    }
end

local H 9
foreach ds in xtarsim_s8y0 xtarsim_s0y8 xtarsim_s8y8 fx_cdols fx_unbalanced {
    local dir = cond(substr("`ds'", 1, 3) == "fx_", "synthetic", "stata")
    import delimited using "$PLP_ROOT/test/fixtures/`dir'/`ds'.csv", clear asdouble case(preserve)
    local bal = "`ds'" != "fx_unbalanced"
    if !`bal' {
        * the unit-level shock x is the CD-OLS shock
        drop s
        rename x s
    }
    xtset id t
    forvalues h = 0/`H' {
        quietly generate double ld`h' = F`h'.y - L.y
    }
    forvalues h = 0/`H' {
        local Y "ld`h'"
        * --- CD-OLS: pooled OLS with an intercept
        quietly regress `Y' s L.s, vce(cluster id)
        _emit `fh' `ds' cdols_cl_id `h'
        quietly regress `Y' s L.s, vce(robust)
        _emit `fh' `ds' cdols_robust `h'
        quietly regress `Y' s L.s LD.y, vce(cluster id)
        _emit `fh' `ds' cdols_dy1_cl_id `h'
        quietly regress `Y' s L.s L2.s, vce(cluster id)
        _emit `fh' `ds' cdols_s2_cl_id `h'
        quietly regress `Y' s L.s i.t, vce(cluster id)
        _emit `fh' `ds' cdols_tfe_cl_id `h'
        quietly regress `Y' s L.s i.t, vce(robust)
        _emit `fh' `ds' cdols_tfe_robust `h'
        if `bal' {
            quietly xtscc `Y' s L.s, lag(3)
            _emit `fh' `ds' cdols_dk3 `h'
        }
        * --- comparator: fixed effects on the cumulative difference (CMFE1)
        quietly xtreg `Y' s L.s, fe vce(cluster id)
        _emit `fh' `ds' cmfe1_cl_id `h'
    }
    * --- locproj, as in the paper: tr(cmlt) builds y_t+h − y_t−1 and Δy lags
    quietly locproj y s, met(reg) tr(cmlt) sl(1) hor(`H') vce(cluster id)
    _emit_irf `fh' `ds' locproj_cmlt_cl_id `H'
    quietly locproj y s, met(reg) tr(cmlt) sl(1) yl(1) hor(`H') vce(cluster id)
    _emit_irf `fh' `ds' locproj_cmlt_yl1_cl_id `H'
    quietly locproj y s, met(reg) tr(cmlt) sl(1) hor(`H') vce(robust)
    _emit_irf `fh' `ds' locproj_cmlt_robust `H'
    * --- comparator: the split-panel jackknife on the level response with xtlp (the SPJ
    * of the paper's applications; its Monte Carlo used locproj m(xtspj), which averages
    * two splits at odd per-unit sample sizes and has a different variance)
    if `bal' {
        quietly xtlp y s, fe m(spj) hor(0 `H') vce(cluster id)
        _emit_hor `fh' `ds' xtlp_spj_level_cl_id `H'
    }
}
file close `fh'
display as result "p3_cdols: done"
