* Phase 1 cross-check: fixed-effects panel LP and its variances.
* Writes test/fixtures/stata/p1_fe.csv (one row per dataset × spec × horizon × term).
* Each spec label maps to a panellp/vcov call in crosscheck/README.md and
* test/test_crosscheck_p1.jl, which rebuilds the derived panels at the end by the same
* filters.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p1_fe.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df_r,N,G,df_a,df_a_initial,df_a_redundant" _n

capture program drop _emit
program define _emit
    * Coefficients and standard errors of the last estimation (e(b), e(V)).
    args fh ds spec h
    foreach v in x L.x L.y w {
        capture scalar __b = _b[`v']
        if _rc continue
        scalar __se = _se[`v']
        file write `fh' "`ds',`spec',`h',`v'," %25.17e (__b) "," %25.17e (__se) "," ///
            %12.0f (e(df_r)) "," %12.0f (e(N)) "," %12.0f (e(N_clust)) "," ///
            %12.0f (e(df_a)) "," %12.0f (e(df_a_initial)) "," %12.0f (e(df_a_redundant)) _n
    }
end

capture program drop _emit_hor
program define _emit_hor
    * Per-horizon results of an xtlp hor(0 H) run (e(b#), e(V#), e(df_r#)).
    args fh ds spec H
    forvalues h = 0/`H' {
        matrix __b = e(b`h')
        matrix __V = e(V`h')
        scalar __dfr = e(df_r`h')
        foreach v in x L.x L.y w {
            local j = colnumb(__b, "`v'")
            if missing(`j') continue
            file write `fh' "`ds',`spec',`h',`v'," %25.17e (__b[1, `j']) "," ///
                %25.17e (sqrt(__V[`j', `j'])) "," %12.0f (__dfr) ",.,.,.,.,." _n
        }
    }
end

capture program drop _leads
program define _leads
    args H
    xtset id t
    forvalues h = 0/`H' {
        quietly generate double F`h'y = F`h'.y
    }
end

local H 6
local X "x L.x L.y w"
foreach ds in fx_balanced fx_unbalanced {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble
    _leads `H'
    forvalues h = 0/`H' {
        local Y "F`h'y"
        * --- xtlp, unit FE
        quietly xtlp `Y' `X', fe m(fe) vce(cluster id)
        _emit `fh' `ds' xtlp_fe_cl_id `h'
        quietly xtlp `Y' `X', fe m(fe) vce(cluster id, ase)
        _emit `fh' `ds' xtlp_fe_cl_id_ase `h'
        quietly xtlp `Y' `X', fe m(fe) vce(robust)
        _emit `fh' `ds' xtlp_fe_robust `h'
        quietly xtlp `Y' `X', fe m(fe) vce(cluster t)
        _emit `fh' `ds' xtlp_fe_cl_t `h'
        quietly xtlp `Y' `X', fe m(fe) vce(cluster id t)
        _emit `fh' `ds' xtlp_fe_cl_idt `h'
        quietly xtlp `Y' `X', fe m(fe) vce(cluster id t, ase)
        _emit `fh' `ds' xtlp_fe_cl_idt_ase `h'
        quietly xtlp `Y' `X', fe m(fe) vce(dkraay)
        _emit `fh' `ds' xtlp_fe_dk `h'
        quietly xtlp `Y' `X', fe m(fe) vce(dkraay lag(3))
        _emit `fh' `ds' xtlp_fe_dk3 `h'
        quietly xtlp `Y' `X', fe m(fe) vce(dkraay lag(3), nodfadj)
        _emit `fh' `ds' xtlp_fe_dk3_nodfadj `h'
        quietly xtlp `Y' `X', fe m(fe) vce(dkraay lag(3), ase)
        _emit `fh' `ds' xtlp_fe_dk3_ase `h'
        quietly xtlp `Y' `X', fe m(fe) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_fe_cl_id_keepsing `h'
        * --- xtlp, unit and time FE
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster id)
        _emit `fh' `ds' xtlp_tfe_cl_id `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster id, ase)
        _emit `fh' `ds' xtlp_tfe_cl_id_ase `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster t)
        _emit `fh' `ds' xtlp_tfe_cl_t `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster id t)
        _emit `fh' `ds' xtlp_tfe_cl_idt `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster id t, ase)
        _emit `fh' `ds' xtlp_tfe_cl_idt_ase `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(robust)
        _emit `fh' `ds' xtlp_tfe_robust `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(dkraay)
        _emit `fh' `ds' xtlp_tfe_dk `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(dkraay lag(3))
        _emit `fh' `ds' xtlp_tfe_dk3 `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(dkraay lag(3), nodfadj)
        _emit `fh' `ds' xtlp_tfe_dk3_nodfadj `h'
        quietly xtlp `Y' `X', tfe m(fe) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_tfe_cl_id_keepsing `h'
        * --- reghdfe (tight tolerance so that iterated demeaning is exact)
        quietly reghdfe `Y' `X', absorb(id t) vce(cluster id t) tolerance(1e-14)
        _emit `fh' `ds' reghdfe_idt_cl_idt `h'
        quietly reghdfe `Y' `X', absorb(id) vce(cluster reg)
        _emit `fh' `ds' reghdfe_id_cl_reg `h'
        quietly reghdfe `Y' `X', absorb(id t grp#t) vce(cluster id) tolerance(1e-14)
        _emit `fh' `ds' reghdfe_3fe_cl_id `h'
        quietly reghdfe `Y' `X', absorb(id) vce(dkraay 4)
        _emit `fh' `ds' reghdfe_id_dk4 `h'
        quietly reghdfe `Y' `X', absorb(id t) vce(dkraay 4) tolerance(1e-14)
        _emit `fh' `ds' reghdfe_idt_dk4 `h'
        quietly reghdfe `Y' `X', absorb(id t) vce(robust) tolerance(1e-14)
        _emit `fh' `ds' reghdfe_idt_robust `h'
        * two-way region × period clusters, not PSD on these panels: the CGM repair fires,
        * with unit effects and in a pooled regression (where _cons is reported)
        quietly reghdfe `Y' `X', absorb(id) vce(cluster reg t)
        _emit `fh' `ds' reghdfe_id_cl_regt `h'
        quietly reghdfe `Y' `X', noabsorb vce(cluster reg t)
        _emit `fh' `ds' reghdfe_noabs_cl_regt `h'
        * --- xtreg (keeps singletons)
        quietly xtreg `Y' `X', fe vce(cluster id)
        _emit `fh' `ds' xtreg_fe_cl_id `h'
        quietly xtreg `Y' `X' i.t, fe vce(cluster id)
        _emit `fh' `ds' xtreg_tfe_cl_id `h'
        * --- pooled OLS
        quietly regress `Y' `X', vce(cluster id)
        _emit `fh' `ds' regress_cl_id `h'
        * --- xtscc refuses data with a period absent for every unit (fx_unbalanced)
        if "`ds'" == "fx_balanced" {
            quietly xtscc `Y' `X', fe lag(3)
            _emit `fh' `ds' xtscc_fe3 `h'
        }
        * --- no lags: the units with one or two rows (fx_unbalanced: 38, 39, 40) enter, and
        * dropping singletons cascades across the unit and time effects
        quietly xtlp `Y' x w, fe m(fe) vce(cluster id)
        _emit `fh' `ds' xtlp_fe_nolag_cl_id `h'
        quietly xtlp `Y' x w, fe m(fe) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_fe_nolag_cl_id_keepsing `h'
        quietly xtlp `Y' x w, tfe m(fe) vce(cluster id)
        _emit `fh' `ds' xtlp_tfe_nolag_cl_id `h'
        quietly reghdfe `Y' x w, absorb(id t) vce(cluster id) tolerance(1e-14)
        _emit `fh' `ds' reghdfe_idt_nolag_cl_id `h'
    }
    * --- response transforms, built by xtlp from y (hor(0 H))
    foreach tr in level diff cmltdiff cmltsum {
        quietly xtlp y `X', fe m(fe) hor(0 `H') ytransf(`tr') vce(cluster id)
        _emit_hor `fh' `ds' xtlp_fe_tr_`tr' `H'
    }
}

* --- Derived panels (rebuilt by the same filters in test/test_crosscheck_p1.jl)

* (a) t ≤ 29: the default Driscoll–Kraay lags ⌊4(T/100)^(2/9)⌋ are 3 at T = 28 (h = 0)
*     and 2 from T = 27 (h ≥ 1), so the T of each horizon's sample is tested
foreach ds in fx_balanced fx_unbalanced {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble
    quietly keep if t <= 29
    _leads `H'
    forvalues h = 0/`H' {
        quietly xtlp F`h'y `X', fe m(fe) vce(dkraay)
        _emit `fh' `ds'_t29 xtlp_fe_dk `h'
        quietly xtlp F`h'y `X', tfe m(fe) vce(dkraay)
        _emit `fh' `ds'_t29 xtlp_tfe_dk `h'
    }
}

* (b) fx_balanced with w missing at t = 20 (a period in the data but in no estimation
*     sample, so calendar and positional lags differ) and unit 30 observed at t = 10, 11
*     only (one complete row at h = 0: a singleton that xtscc keeps)
import delimited using "$PLP_ROOT/test/fixtures/synthetic/fx_balanced.csv", clear asdouble
quietly drop if id == 30 & !inlist(t, 10, 11)
quietly replace w = . if t == 20
_leads `H'
forvalues h = 0/`H' {
    quietly xtscc F`h'y `X', fe lag(3)
    _emit `fh' fx_balanced_gap xtscc_fe3 `h'
    quietly xtlp F`h'y `X', fe m(fe) vce(dkraay lag(3))
    _emit `fh' fx_balanced_gap xtlp_fe_dk3 `h'
}

* (c) two disconnected blocks: units 1–15 at t ≤ 30 and units 16–30 at t > 30, so the
*     unit and time effects have two connected components
import delimited using "$PLP_ROOT/test/fixtures/synthetic/fx_balanced.csv", clear asdouble
quietly drop if (id <= 15 & t > 30) | (id > 15 & t <= 30)
_leads `H'
forvalues h = 0/`H' {
    quietly reghdfe F`h'y `X', absorb(id t) vce(cluster id) tolerance(1e-14)
    _emit `fh' fx_balanced_2blocks reghdfe_idt_cl_id `h'
    quietly reghdfe F`h'y `X', absorb(id t) vce(robust) tolerance(1e-14)
    _emit `fh' fx_balanced_2blocks reghdfe_idt_robust `h'
}

file close `fh'
display as result "p1_fe: done"
