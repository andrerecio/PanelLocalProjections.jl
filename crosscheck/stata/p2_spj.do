* Phase 2 cross-check: split-panel jackknife (Mei, Sheng & Shi 2026) with xtlp m(spj), and
* the xtspj point estimates (Sun & Dhaene 2019) on balanced panels with an even number of
* periods. Writes test/fixtures/stata/p2_spj.csv.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p2_spj.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df_r,N,G,df_a" _n

capture program drop _emit
program define _emit
    args fh ds spec h
    foreach v in x L.x L.y w Lx Ly {
        capture scalar __b = _b[`v']
        if _rc continue
        scalar __se = _se[`v']
        local name = cond("`v'" == "Lx", "L.x", cond("`v'" == "Ly", "L.y", "`v'"))
        file write `fh' "`ds',`spec',`h',`name'," %25.17e (__b) "," %25.17e (__se) "," ///
            %12.0f (e(df_r)) "," %12.0f (e(N)) "," %12.0f (e(N_clust)) "," %12.0f (e(df_a)) _n
    }
end

foreach ds in fx_balanced fx_unbalanced {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble
    xtset id t
    forvalues h = 0/3 {
        quietly generate double F`h'y = F`h'.y
    }
    local X "x L.x L.y w"
    forvalues h = 0/3 {
        local Y "F`h'y"
        quietly xtlp `Y' `X', fe m(spj) vce(cluster id, ase)
        _emit `fh' `ds' xtlp_spj_fe_cl_id_ase `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster id)
        _emit `fh' `ds' xtlp_spj_fe_cl_id `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster id t, ase)
        _emit `fh' `ds' xtlp_spj_fe_cl_idt_ase `h'
        quietly xtlp `Y' `X', fe m(spj) vce(dkraay lag(2), ase)
        _emit `fh' `ds' xtlp_spj_fe_dk2_ase `h'
        quietly xtlp `Y' `X', fe m(spj) vce(robust)
        _emit `fh' `ds' xtlp_spj_fe_robust `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster id, ase) keepsingletons
        _emit `fh' `ds' xtlp_spj_fe_cl_id_ase_keepsing `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster t)
        _emit `fh' `ds' xtlp_spj_fe_cl_t `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster id t)
        _emit `fh' `ds' xtlp_spj_fe_cl_idt `h'
        quietly xtlp `Y' `X', fe m(spj) vce(dkraay)
        _emit `fh' `ds' xtlp_spj_fe_dk `h'
        quietly xtlp `Y' `X', fe m(spj) vce(dkraay lag(2), nodfadj)
        _emit `fh' `ds' xtlp_spj_fe_dk2_nodfadj `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(cluster id, ase)
        _emit `fh' `ds' xtlp_spj_tfe_cl_id_ase `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(cluster id)
        _emit `fh' `ds' xtlp_spj_tfe_cl_id `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(cluster t)
        _emit `fh' `ds' xtlp_spj_tfe_cl_t `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(cluster id t, ase)
        _emit `fh' `ds' xtlp_spj_tfe_cl_idt_ase `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(robust)
        _emit `fh' `ds' xtlp_spj_tfe_robust `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(dkraay lag(2), ase)
        _emit `fh' `ds' xtlp_spj_tfe_dk2_ase `h'
        quietly xtlp `Y' `X', tfe m(spj) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_spj_tfe_cl_id_keepsing `h'
        quietly xtlp `Y' `X', fe m(spj) vce(cluster reg)
        _emit `fh' `ds' xtlp_spj_fe_cl_reg `h'
        * no lags: units with two or three rows enter, with single-row halves
        quietly xtlp `Y' x w, fe m(spj) vce(cluster id)
        _emit `fh' `ds' xtlp_spj_fe_nolag_cl_id `h'
        quietly xtlp `Y' x w, fe m(spj) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_spj_fe_nolag_cl_id_keepsing `h'
        quietly xtlp `Y' x w, tfe m(spj) vce(cluster id)
        _emit `fh' `ds' xtlp_spj_tfe_nolag_cl_id `h'
        quietly xtlp `Y' x w, tfe m(spj) vce(cluster id) keepsingletons
        _emit `fh' `ds' xtlp_spj_tfe_nolag_cl_id_keepsing `h'
    }
    * xtspj needs every unit to have the same even number of contiguous periods: on the
    * balanced panel with one lag, h = 1 and h = 3 leave 58 and 56 periods per unit.
    if "`ds'" == "fx_balanced" {
        quietly generate double Lx = L.x
        quietly generate double Ly = L.y
        foreach h in 1 3 {
            preserve
            quietly keep if !missing(F`h'y, x, Lx, Ly, w)
            quietly bysort id: assert _N == 59 - `h'
            quietly bysort id: assert mod(_N, 2) == 0
            quietly xtspj F`h'y x Lx Ly w, model(regress) method(parm)
            _emit `fh' `ds' xtspj_parm `h'
            restore
        }
    }
}
file close `fh'
display as result "p2_spj: done"
