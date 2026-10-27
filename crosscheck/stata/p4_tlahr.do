* Phase 4 cross-check: t-LAHR inference (Almuzara & Sancibrián 2026) on the fx_tlahr panels.
*
* A. Lag-augmented two-way FE regressions (reghdfe): y_i,t+h on s_i X_t, p = min(h, 2) lags
*    of y and of s X, unit and time effects, clustered by time. Both panels.
* B. The synthetic regression (supplement B.2) on the balanced panel, where the panel
*    regression *is* a time-series regression: every regressor is s_i times an aggregate
*    series (s ≡ 1 with unit effects; s_i with unit and time effects), so the t-LAHR
*    variance with the HC2/Imbens–Kolesár refinement is `regress …, noconstant
*    vce(hc2, dfadjust)` on the ŝ-weighted cross-sectional averages, each column demeaned
*    over the sample periods, and the equal-weighted cosine variance is `harreg`.
* Writes test/fixtures/stata/p4_tlahr.csv; spec labels map to panellp/vcov calls in
* test/test_crosscheck_p4.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p4_tlahr.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df,N,G" _n

local H 6
* --- A. lag-augmented regressions, CR1 clustered by time (reghdfe)
foreach ds in fx_tlahr fx_tlahr_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    generate double sX = s * X
    forvalues h = 0/`H' {
        quietly generate double Fy`h' = F`h'.y
        local p = min(`h', 2)
        local lags ""
        if `p' > 0 local lags "L(1/`p').y L(1/`p').sX"
        quietly reghdfe Fy`h' sX `lags', absorb(id t) vce(cluster t)
        file write `fh' "`ds',la_cap2_twfe_cl_t,`h',sX," %25.17e (_b[sX]) "," ///
            %25.17e (_se[sX]) "," %12.0f (e(df_r)) "," %12.0f (e(N)) "," ///
            %12.0f (e(N_clust)) _n
    }
}

* --- B. synthetic regressions on the balanced panel, p = 2 lags of the shock
foreach spec in ones_unit s_twfe {
    forvalues h = 0/`H' {
        import delimited using "$PLP_ROOT/test/fixtures/synthetic/fx_tlahr.csv", clear asdouble case(preserve)
        xtset id t
        generate double Fy = F`h'.y
        * (sum) of an all-missing period would be 0: drop the rows without the lead
        quietly drop if missing(Fy)
        * weights: 1 (s ≡ 1, unit effects) or ŝ_i = s_i − mean(s) (time effects)
        if "`spec'" == "ones_unit" generate double w = 1
        else {
            quietly summarize s
            generate double w = s - r(mean)
        }
        generate double wy = w * Fy
        generate double w2 = w * w
        collapse (sum) wy w2 (mean) X, by(t)
        generate double yw = wy / w2
        tsset t
        generate double L1X = L.X
        generate double L2X = L2.X
        * the sample: periods with the lead and both lags
        quietly keep if !missing(yw, X, L1X, L2X)
        foreach v in yw X L1X L2X {
            quietly summarize `v'
            quietly replace `v' = `v' - r(mean)
        }
        quietly regress yw X L1X L2X, noconstant vce(hc2, dfadjust)
        matrix __df = e(adj_df)
        file write `fh' "fx_tlahr,syn_`spec'_hc2,`h',X," %25.17e (_b[X]) "," ///
            %25.17e (_se[X]) "," %25.17e (__df[1, colnumb(__df, "X")]) "," ///
            %12.0f (e(N)) ",." _n
        quietly harreg yw X L1X L2X, noconstant estimator(ewc)
        file write `fh' "fx_tlahr,syn_`spec'_ewc,`h',X," %25.17e (_b[X]) "," ///
            %25.17e (_se[X]) "," %25.17e (e(df_fb)) "," %12.0f (e(N)) ",." _n
    }
}
file close `fh'
display as result "p4_tlahr: done"
