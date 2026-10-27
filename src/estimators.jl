# Per-horizon fits, dispatched on the estimator type.

# A fit without coefficients: no observations (:empty), or a shock that is not identified
# in the horizon's sample (:singular, :spj_half_singular), with `nobs` rows in it.
function _empty_fit(d::PanelDesign, h::Int, mask::BitVector, singletons::Int;
        status::Symbol = :empty, nobs::Int = 0, keep::BitVector = falses(length(d.coefnames)),
        absorbed::BitVector = falses(length(d.coefnames)), iterations::Int = 0,
        converged::Bool = true, active::BitVector = trues(length(d.coefnames)),
        nlags::Int = 0, instkeep::BitVector = falses(size(d.Z, 2)))
    K = length(d.coefnames)
    G = length(d.groups)
    D = length(d.fe)
    dof = DofInfo(nobs, 0, 0, zeros(Int, D), zeros(Int, D, D), zeros(Int, G),
        Matrix{Bool}(falses(D, G)), singletons)
    stats = ScoreStats(zeros(K, K), [zeros(ng, K) for ng in d.ngroups],
        Tuple{Int, Int}[], Matrix{Float64}[], Symbol[])
    return HorizonFit(h, status, fill(NaN, K), keep, absorbed, fill(NaN, K, K),
        fill(NaN, K), mask, falses(d.ngroups[2]), dof, stats, iterations, converged,
        active, nlags, SyntheticSeries[], instkeep, nothing)
end

"""
    _within_fit(d, y, mask, absorb, cols = all columns, zcols = all instruments) -> NamedTuple

Within transformation and OLS on the rows of `mask` and the columns `cols` (ascending, the
shocks among them). `keep`, `absorbed` and `xsd` are over all columns of the design;
`Xk` holds the within-transformed kept columns, `β` and `B` their coefficients and bread.
The design is absorbed in place in one buffer, of which `Xk` is a view: the raw columns
are not kept (their statistics are taken before the absorption).

With an IV block (`d.endo`), the excluded instruments are absorbed in the same buffer and
the fit is two-stage least squares: the first stage projects the endogenous columns on the
exogenous ones and the instruments (`_firststage_basis`), whose fitted values X̂ replace
them in `Xk`; `B` is (X̂'X̂)⁻¹ and `resid` the structural residual ỹ − X̃β (with the
endogenous columns themselves), so that the scores are X̂ ⊙ resid. `zcols` are the
instruments of this fit (`LagAugmented`: those with lags ≤ p(h)), and `instkeep` (over all
instruments) marks those in the first stage; `status` is `:underidentified` when fewer of them than
endogenous columns are left.
"""
function _within_fit(d::PanelDesign, y::Vector{Float64}, mask::BitVector,
        absorb::AbsorptionMethod, cols::Vector{Int} = collect(axes(d.X, 2)),
        zcols::Vector{Int} = collect(axes(d.Z, 2)))
    rows = findall(mask)
    nc, L = length(cols), length(zcols)
    A = Matrix{Float64}(undef, length(rows), 1 + nc + L)
    A[:, 1] = view(y, rows)
    for (c, j) in enumerate(cols)
        A[:, c + 1] = view(d.X, rows, j)
    end
    for (l, k) in enumerate(zcols)
        A[:, 1 + nc + l] = view(d.Z, rows, k)
    end
    rawall = _colstats(view(A, :, 2:size(A, 2)))
    raw = (; rawall.n, css = rawall.css[1:nc], ss = rawall.ss[1:nc])
    fs = _fe_sample(d.fe, rows)
    it, conv = _absorb_matrix!(A, fs, absorb)
    yt = A[:, 1]
    hasfe = !isempty(d.fe)
    shocks = Int[something(findfirst(==(j), cols)) for j in d.shockcols]
    # lags added by LagAugmented give way to the formula's columns when collinear
    added = Int[k for (k, j) in enumerate(cols) if d.lagorder[j] > 0]
    Xt = view(A, :, 2:(1 + nc))
    absall = _absorbed_cols(view(A, :, 2:size(A, 2)), rawall, hasfe, d.intercept)
    xabs = absall[1:nc]
    K = length(d.coefnames)
    instkeep = falses(size(d.Z, 2))
    instkeep[zcols] = .!absall[(nc + 1):end]
    # 2SLS: the fitted values of the (kept) endogenous columns replace them in the buffer
    en = Int[k for (k, j) in enumerate(cols) if j in d.endo && !xabs[k]]
    Xen = zeros(length(rows), 0)
    fst = (; endo = Int[], inst = Int[], coef = zeros(0, 0), zz = zeros(0, 0),
        yy = zeros(0, 0), U = zeros(0, 0))
    if L > 0 && !isempty(en)
        ex = Int[k for k in 1:nc if !xabs[k] && !(k in en)]
        zi = findall(.!absall[(nc + 1):end])
        Q, zk, rE = _firststage_basis(view(A, :, 2:size(A, 2)), ex, nc .+ zi)
        instkeep[zcols[zi]] = zk
        Xen = A[:, 1 .+ en]
        if count(zk) < length(en) || !_identified(Q, rE, Xen)
            keep = falses(K)
            absorbed = falses(K)
            absorbed[cols] = xabs
            return (; β = Float64[], B = zeros(0, 0), keep, absorbed, resid = Float64[],
                status = :underidentified, rows, Xk = view(A, :, 2:1), yt, fs,
                xsd = fill(NaN, K), iterations = it, converged = conv, instkeep, first = fst)
        end
        fsr = _firststage_stats(view(A, :, 2:size(A, 2)), Q, rE, Xen, nc .+ zi[zk])
        fst = (; endo = cols[en], inst = zcols[zi[zk]], fsr.coef, fsr.zz, fsr.yy, fsr.U)
        A[:, 1 .+ en] = Q * (Q' * Xen)
    end
    o = _ols(Xt, yt, raw, hasfe, d.intercept, shocks, added; absorbed = xabs)
    resid = o.resid
    if o.status === :ok && !isempty(en)
        # ỹ − X̃β = (ỹ − X̂β) − (X̃ − X̂)β, the difference only in the endogenous columns
        kidx = findall(o.keep)
        for (c, k) in enumerate(en)
            o.keep[k] || continue
            b = o.β[findfirst(==(k), kidx)]
            resid .-= b .* (view(Xen, :, c) .- view(A, :, 1 + k))
        end
    end
    keep, absorbed, xsd = falses(K), falses(K), fill(NaN, K)
    keep[cols] = o.keep
    absorbed[cols] = o.absorbed
    xsd[cols] = _xsd(raw, d.intercept)
    # the kept columns, moved to the front of the buffer (a contiguous view)
    kc = findall(o.keep)
    for (c, k) in enumerate(kc)
        c == k || (A[:, c + 1] = view(A, :, k + 1))
    end
    return (; o.β, o.B, keep, absorbed, resid, o.status, rows,
        Xk = view(A, :, 2:(1 + length(kc))), yt, fs, xsd, iterations = it,
        converged = conv, instkeep, first = fst)
end

"""
    _firststage_stats(M, Q, rE, Xen, zc) -> NamedTuple

The first stage of 2SLS for the diagnostics of `firststage`: the endogenous columns `Xen`
and the instruments (columns `zc` of `M`), both net of the exogenous columns (the first
`rE` columns of the basis `Q`), the coefficients Π̂ of the regressions of the former on the
latter, the cross-products, and the scores z̃ v̂_k (k outer).
"""
function _firststage_stats(M::AbstractMatrix{Float64}, Q::Matrix{Float64}, rE::Int,
        Xen::Matrix{Float64}, zc::Vector{Int})
    QE = view(Q, :, 1:rE)
    Z = M[:, zc]
    Z .-= QE * (QE' * Z)
    Y = Xen .- QE * (QE' * Xen)
    zz = Z' * Z
    coef = zz \ (Z' * Y)
    V = Y .- Z * coef
    L = size(Z, 2)
    U = Matrix{Float64}(undef, size(Z, 1), L * size(Y, 2))
    for k in axes(Y, 2)
        view(U, :, ((k - 1) * L + 1):(k * L)) .= Z .* view(V, :, k)
    end
    return (; coef, zz, yy = Y' * Y, U)
end

"""
    _firststage_basis(M, ex, zi) -> (Q, zkeep, rE)

First stage of 2SLS on the within-transformed columns of `M`: an orthonormal basis `Q` of
the span of the exogenous columns `ex` and the instruments `zi`, and which instruments
enter it. Frisch–Waugh–Lovell, as `_dropcollinear!`: the pivoted QR of the exogenous
columns (unit norm, rank at the tolerance of `_ols`); then each instrument in turn (unit
norm), residualised on the basis so far (twice, for orthogonality to round-off), enters
when the residual norm exceeds the same tolerance. An instrument collinear with the
exogenous columns or with earlier instruments is left out (as ivreg2 drops the later of
collinear instruments). The first `rE` columns of `Q` span the exogenous columns.
"""
function _firststage_basis(M::AbstractMatrix{Float64}, ex::Vector{Int}, zi::Vector{Int})
    n = size(M, 1)
    QE = zeros(n, 0)
    if !isempty(ex)
        F = LinearAlgebra.qr!(_unitcols!(M[:, ex]), LinearAlgebra.ColumnNorm())
        dg = abs.(LinearAlgebra.diag(F.R))
        r = count(>(_RANK_TOL * dg[1]), dg)
        QE = F.Q * Matrix{Float64}(LinearAlgebra.I, n, r)
    end
    R = _unitcols!(M[:, zi])
    zkeep = falses(length(zi))
    QZ = Vector{Float64}[]
    for c in eachindex(zi)
        v = R[:, c]
        for _ in 1:2
            v .-= QE * (QE' * v)
            for q in QZ
                v .-= LinearAlgebra.dot(q, v) .* q
            end
        end
        nv = LinearAlgebra.norm(v)
        nv > _RANK_TOL || continue
        zkeep[c] = true
        push!(QZ, v ./ nv)
    end
    return isempty(QZ) ? QE : hcat(QE, reduce(hcat, QZ)), zkeep, size(QE, 2)
end

"""
    _identified(Q, rE, X) -> Bool

Whether the instruments identify the endogenous columns `X` (within-transformed): the
parts of the columns off the exogenous span (the first `rE` columns of the first-stage
basis `Q`) that the instruments (its other columns) predict have full column rank. A
column whose part predicted by the instruments is round-off, a partial R² of the first
stage at most `_ABSORBED_TOL` (the absorbed-column screen's tolerance), is not
identified; nor are columns whose predicted parts are collinear. Its fitted value would
otherwise be round-off rescaled into an arbitrary coefficient. Columns spanned by the
exogenous regressors are left to the collinearity screen of `_ols`.
"""
function _identified(Q::Matrix{Float64}, rE::Int, X::Matrix{Float64})
    QE = view(Q, :, 1:rE)
    QZ = view(Q, :, (rE + 1):size(Q, 2))
    R = X .- QE * (QE' * X)
    cols = [c
            for c in axes(X, 2)
            if sum(abs2, view(R, :, c)) > _ABSORBED_TOL * sum(abs2, view(X, :, c))]
    isempty(cols) && return true
    size(QZ, 2) >= length(cols) || return false
    G = QZ' * _unitcols!(R[:, cols])
    return minimum(LinearAlgebra.svdvals(G)) > sqrt(_ABSORBED_TOL)
end

# The columns of C scaled to unit norm (zero columns left as they are).
function _unitcols!(C::Matrix{Float64})
    for c in axes(C, 2)
        nc = LinearAlgebra.norm(view(C, :, c))
        nc > 0 && (C[:, c] ./= nc)
    end
    return C
end

function _singular_fit(d::PanelDesign, h::Int, mask::BitVector, singletons::Int, w;
        status::Symbol = :singular, active::BitVector = trues(length(d.coefnames)),
        nlags::Int = 0)
    return _empty_fit(d, h, mask, singletons; status, nobs = length(w.rows),
        keep = w.keep, absorbed = w.absorbed, iterations = w.iterations,
        converged = w.converged, active, nlags, instkeep = w.instkeep)
end

function _dofinfo(d::PanelDesign, rows::Vector{Int}, fs::_FESample, k::Int,
        present::Vector{Int}, singletons::Int)
    L = length.(fs.counts)
    m = _fe_components(fs)
    df_a = _df_fe(L, m, collect(eachindex(L)), false)
    nested = Matrix{Bool}(undef, length(d.fe), length(d.groups))
    for f in eachindex(d.fe), g in eachindex(d.groups)

        nested[f, g] = _isnested(fs.codes[f], L[f], view(d.groups[g], rows))
    end
    return DofInfo(length(rows), k, df_a, L, m, present, nested, singletons)
end

# Standard deviations (n − 1) of the raw regressors, floored at 1e-3: the metric in which
# reghdfe and xtlp clip negative eigenvalues of a multi-way variance. The intercept enters
# with scale 1, as `_cons` in reghdfe (reghdfe.mata, `stdev_x = stdev_x, 1`).
function _xsd(raw::NamedTuple, intercept::Bool)
    n = raw.n
    sd = [n > 1 ? max(sqrt(css / (n - 1)), 1e-3) : 1e-3 for css in raw.css]
    intercept && (sd[1] = 1.0)
    return sd
end

function _expand(β::Vector{Float64}, B::Matrix{Float64}, keep::BitVector)
    K = length(keep)
    kidx = findall(keep)
    coef = fill(NaN, K)
    coef[kidx] = β
    bread = fill(NaN, K, K)
    bread[kidx, kidx] = B
    return coef, bread
end

"""
    _fit_horizon(estimator, design, h, y, mask, singletons, absorb) -> HorizonFit
"""
function _fit_horizon(::WithinFE, d::PanelDesign, h::Int, y::Vector{Float64},
        mask::BitVector, singletons::Int, absorb::AbsorptionMethod)
    return _fit_within(d, h, y, mask, singletons, absorb, collect(axes(d.X, 2)), 0)
end

# Fixed-effects fit of horizon h on the columns `cols`, with the synthetic series of the
# shocks.
function _fit_within(d::PanelDesign, h::Int, y::Vector{Float64}, mask::BitVector,
        singletons::Int, absorb::AbsorptionMethod, cols::Vector{Int}, nlags::Int)
    active = falses(length(d.coefnames))
    active[cols] .= true
    any(mask) || return _empty_fit(d, h, mask, singletons; active, nlags)
    w = _within_fit(d, y, mask, absorb, cols, findall(<=(nlags), d.zlagorder))
    w.status === :ok ||
        return _singular_fit(d, h, mask, singletons, w; status = w.status, active, nlags)
    groups = [view(g, w.rows) for g in d.groups]
    timepresent = _present(groups[2], d.ngroups[2])
    # IV fits have no synthetic regression (TLAHR(refine = true) is not defined for them)
    synthetic = isempty(d.endo) ? _synthetic(d, w, groups[2], timepresent) :
                SyntheticSeries[]
    # the scores, in the buffer of the within-transformed design (no longer needed)
    U = w.Xk
    U .*= w.resid
    stats, present = _score_stats(U, w.keep, groups, d.ngroups)
    dof = _dofinfo(d, w.rows, w.fs, count(w.keep), present, singletons)
    coef, bread = _expand(w.β, w.B, w.keep)
    fstage = nothing
    if !isempty(w.first.endo)
        m = length(w.first.endo)
        fstats, _ = _score_stats(w.first.U, trues(size(w.first.U, 2)), groups, d.ngroups)
        fstage = FirstStage(w.first.endo, w.first.inst, w.first.coef,
            kron(Matrix{Float64}(LinearAlgebra.I, m, m), inv(w.first.zz)), fstats,
            w.first.zz, w.first.yy)
    end
    return HorizonFit(h, :ok, coef, w.keep, w.absorbed, bread, w.xsd, mask, timepresent,
        dof, stats, w.iterations, w.converged, active, nlags, synthetic, w.instkeep, fstage)
end

"""
    _shockweights(d, j, rows, keep) -> (w, raw)

Weights of the synthetic regression of column j (an aggregate series times an exposure)
on the rows `rows`: its exposure residualised per period on an intercept when a fixed
effect nested in the periods is absorbed (demeaned within its groups: Almuzara &
Sancibrián's ŝ_i), and on the exposures of the other kept regressors (`keep`) that
multiply the same aggregate series, shocks or not (by Frisch–Waugh–Lovell, the variation
that identifies coefficient j); and the raw exposure. A period whose residualised exposure
is round-off (sum of squares at most 1e-20 of the raw exposure's, e.g. no more units than
exposures to residualise on) has weight zero.
"""
function _shockweights(d::PanelDesign, j::Int, rows::Vector{Int},
        keep::AbstractVector{Bool})
    e = d.exposures[j][rows]
    others = Int[l
                 for l in eachindex(d.exposures)
                 if l != j && keep[l] && _sameagg(d.aggseries[l], d.aggseries[j])]
    isempty(d.timefe) && isempty(others) && return e, e
    A = Matrix{Float64}(undef, length(rows), 1 + length(others))
    A[:, 1] = e
    for (c, l) in enumerate(others)
        A[:, c + 1] = view(d.exposures[l], rows)
    end
    isempty(d.timefe) ||
        _absorb_matrix!(A, _fe_sample(d.fe[d.timefe], rows), AlternatingProjections())
    w = A[:, 1]
    per = view(d.pidx.period, rows)
    isempty(others) || _perperiod_resid!(w, view(A, :, 2:size(A, 2)), per)
    ss = zeros(nspan(d.pidx), 2)
    for i in eachindex(w, e)
        k = per[i] - d.pidx.pmin + 1
        ss[k, 1] += w[i]^2
        ss[k, 2] += e[i]^2
    end
    for i in eachindex(w)
        k = per[i] - d.pidx.pmin + 1
        ss[k, 1] <= 1e-20 * ss[k, 2] && (w[i] = 0.0)
    end
    return w, e
end

# w minus its projection on the columns of Z, period by period (`per`).
function _perperiod_resid!(w::Vector{Float64}, Z::AbstractMatrix{Float64},
        per::AbstractVector{<:Integer})
    order = sortperm(per)
    n = length(order)
    i = 1
    while i <= n
        k = i
        while k < n && per[order[k + 1]] == per[order[i]]
            k += 1
        end
        rr = order[i:k]
        _projresid!(view(w, rr), Z[rr, :])
        i = k + 1
    end
    return w
end

# w minus its projection on the span of the columns of Z: a pivoted QR of the columns
# scaled to unit norm, of rank at the tolerance of `_ols` (never a solve, which a singular
# block would make fail).
function _projresid!(w::AbstractVector{Float64}, Z::Matrix{Float64})
    for c in axes(Z, 2)
        nc = LinearAlgebra.norm(view(Z, :, c))
        nc > 0 && (Z[:, c] ./= nc)
    end
    F = LinearAlgebra.qr!(Z, LinearAlgebra.ColumnNorm())
    r = count(>(_RANK_TOL), abs.(LinearAlgebra.diag(F.R)))
    r == 0 && return w
    Q = F.Q * Matrix{Float64}(LinearAlgebra.I, size(Z, 1), r)
    w .-= Q * (Q' * w)
    return w
end

# Σ_{r∈t} w_r z_r per period (rows of `tpos`, 1:T) for each column of Z.
function _persums(Z::AbstractVecOrMat{Float64}, w::Vector{Float64}, tpos::Vector{Int},
        T::Int)
    S = zeros(T, size(Z, 2))
    @inbounds for c in axes(Z, 2), r in axes(Z, 1)

        S[tpos[r], c] += w[r] * Z[r, c]
    end
    return S
end

# Per-period weighted means Σ w z / Σ w², over the periods present; `zero` where a period
# has no weight.
function _wmeans(Z, w::Vector{Float64}, tpos::Vector{Int}, ω::Vector{Float64}, zero)
    S = _persums(Z, w, tpos, length(ω))
    for t in eachindex(ω)
        S[t, :] .= ω[t] > 0 ? S[t, :] ./ ω[t] : zero
    end
    return S
end

# Synthetic series of every kept coefficient that owns a synthetic regression
# (`PanelDesign.owner`; Almuzara & Sancibrián 2026, Remark 1; supplement B.2), from a fit's
# within-transformed data. Columns that the weights annihilate (the other exposures of the
# same aggregate series, on balanced panels) are round-off, which the rank tolerance of
# `_ikdesign` drops.
function _synthetic(d::PanelDesign, w, tcodes::AbstractVector{Int32},
        timepresent::BitVector)
    out = SyntheticSeries[]
    pos = zeros(Int, length(timepresent))
    pos[timepresent] = 1:count(timepresent)
    tpos = [pos[c] for c in tcodes]
    T = count(timepresent)
    K = length(d.coefnames)
    kidx = findall(w.keep)
    for j in sort!(unique(d.owner))
        (j > 0 && w.keep[j]) || continue
        wt, raw = _shockweights(d, j, w.rows, w.keep)
        ω = vec(_persums(wt, wt, tpos, T))
        ωraw = vec(_persums(raw, raw, tpos, T))
        X = fill(NaN, T, K)
        X[:, kidx] = _wmeans(w.Xk, wt, tpos, ω, 0.0)
        Xraw = fill(NaN, T, K)
        Xraw[:, kidx] = _wmeans(w.Xk, raw, tpos, ωraw, 0.0)
        push!(out,
            SyntheticSeries(j, ω, X, vec(_wmeans(w.yt, wt, tpos, ω, NaN)),
                vec(_wmeans(w.resid, wt, tpos, ω, NaN)), Xraw))
    end
    return out
end

# Half (1 or 2) of each estimation row: per unit, rows in time order, the first ⌈n_i/2⌉
# (or ⌊n_i/2⌋) in half 1. Mei–Sheng–Shi eq. (16) and appendix B.2.4; LP_panel_all.r and
# xtlp split by count within each unit, which is `:ceil`.
function _spj_halves(pidx::PanelIndex, rows::Vector{Int}, split::Symbol)
    order = sortperm(rows; by = r -> (pidx.unit[r], pidx.period[r]))
    half = Vector{Int8}(undef, length(rows))
    n = length(rows)
    k = 1
    while k <= n
        u = pidx.unit[rows[order[k]]]
        stop = k
        while stop < n && pidx.unit[rows[order[stop + 1]]] == u
            stop += 1
        end
        ni = stop - k + 1
        na = split === :ceil ? cld(ni, 2) : fld(ni, 2)
        for j in k:stop
            half[order[j]] = j - k < na ? Int8(1) : Int8(2)
        end
        k = stop + 1
    end
    return half
end

# OLS on a split-panel half; `nothing` when the half's design is rank-deficient.
function _half_ols(X::Matrix{Float64}, y::Vector{Float64})
    size(X, 1) >= size(X, 2) || return nothing
    norms = [LinearAlgebra.norm(view(X, :, j)) for j in axes(X, 2)]
    any(iszero, norms) && return nothing
    F = LinearAlgebra.qr(X ./ norms', LinearAlgebra.ColumnNorm())
    dg = abs.(LinearAlgebra.diag(F.R))
    count(>(_RANK_TOL * dg[1]), dg) == size(X, 2) || return nothing
    return (F \ y) ./ norms
end

function _fit_horizon(est::SplitPanelJackknife, d::PanelDesign, h::Int,
        y::Vector{Float64}, mask::BitVector, singletons::Int, absorb::AbsorptionMethod)
    any(mask) || return _empty_fit(d, h, mask, singletons)
    any(c -> c == d.pidx.unit, d.fe) || throw(ArgumentError(
        "SplitPanelJackknife needs a unit fixed effect, fe(<id column>), in the formula"))
    w = _within_fit(d, y, mask, absorb)
    w.status === :singular && return _singular_fit(d, h, mask, singletons, w)
    kidx = findall(w.keep)
    half = _spj_halves(d.pidx, w.rows, est.split)
    Xt = w.Xk
    Xhalf = similar(Xt)
    θ = Vector{Vector{Float64}}(undef, 2)
    iterations, converged = w.iterations, w.converged
    for k in 1:2
        sel = findall(==(Int8(k)), half)
        rk = w.rows[sel]
        Xr = d.X[rk, kidx]
        A = hcat(y[rk], Xr)
        it, conv = _absorb_matrix!(A, _fe_sample(d.fe, rk), absorb)
        iterations = max(iterations, it)
        converged &= conv
        # a regressor absorbed by the half's own fixed effects leaves a round-off column
        # that rescaling would turn into an arbitrary coefficient: the half is singular
        ab = _absorbed_cols(view(A, :, 2:size(A, 2)), Xr, true, false)
        θk = any(ab) ? nothing : _half_ols(A[:, 2:end], A[:, 1])
        if θk === nothing
            # the fit's `absorbed` marks the kept columns absorbed inside the half
            habs = falses(length(w.keep))
            habs[kidx[ab]] .= true
            return _singular_fit(d, h, mask, singletons,
                merge(w, (; iterations, converged, absorbed = habs));
                status = :spj_half_singular)
        end
        θ[k] = θk
        Xhalf[sel, :] = view(A, :, 2:size(A, 2))
    end
    # θ̂_spj = 2θ̂ − ½(θ̂_a + θ̂_b); residuals on the full-sample within data; scores
    # d = 2X̃ − X̃_half (LP_panel_all.r:893–923; xtlp.ado:1369–1412)
    βspj = 2 .* w.β .- (θ[1] .+ θ[2]) ./ 2
    resid = w.yt .- Xt * βspj
    U = (2 .* Xt .- Xhalf) .* resid
    groups = [view(g, w.rows) for g in d.groups]
    stats, present = _score_stats(U, w.keep, groups, d.ngroups)
    dof = _dofinfo(d, w.rows, w.fs, count(w.keep), present, singletons)
    coef, bread = _expand(βspj, w.B, w.keep)
    return HorizonFit(h, :ok, coef, w.keep, w.absorbed, bread, w.xsd, mask,
        _present(groups[2], d.ngroups[2]), dof, stats, iterations, converged,
        trues(length(w.keep)), 0, SyntheticSeries[], falses(0), nothing)
end

# ---------------------------------------------------------------------------------------
# Specification checks, run by `panellp` before estimation
# ---------------------------------------------------------------------------------------

_validate(::PanelLPEstimator, d::PanelDesign, f) = nothing

# Panel LP-IV (Almuzara & Sancibrián 2026, §3.4) is two-stage least squares on the
# fixed-effects design, with or without lag augmentation; the other estimators have no IV
# version.
_check_iv(::Union{WithinFE, LagAugmented}) = nothing
function _check_iv(est::PanelLPEstimator)
    why = est isa SplitPanelJackknife ?
          "the split-panel jackknife (Mei, Sheng & Shi 2026) is defined for OLS, and its " *
          "bias correction and variance are not established for 2SLS" :
          est isa CDOLS ?
          "cumulative-difference OLS (Ugarte-Ruiz 2026) has no IV version" :
          "$(nameof(typeof(est))) has no IV version"
    throw(ArgumentError(
        "an instrumental-variable block (endo ~ inst) needs estimator = WithinFE() or " *
        "LagAugmented(…): $why"))
end

# Every use of the variables `ys` in term `t`, as (smallest lag, differenced).
_uses!(out, t::VarTerm, ys, lag, diff) = (t.sym in ys && push!(out, (lag, diff)); out)
_uses!(out, ::ConstTerm, ys, lag, diff) = out
_uses!(out, t::LagsTerm, ys, lag, diff) = _uses!(out, t.term, ys, lag + first(t.lags), diff)
_uses!(out, t::FirstDiffTerm, ys, lag, diff) = _uses!(out, t.term, ys, lag, true)
function _uses!(out, t::FunTerm, ys, lag, diff)
    foreach(a -> _uses!(out, a, ys, lag, diff), t.args)
    return out
end
function _uses!(out, t::ProductTerm, ys, lag, diff)
    foreach(a -> _uses!(out, a, ys, lag, diff), t.terms)
    return out
end

# Date of each column of a term: the lag ℓ of its most recent value, x_{t−ℓ} (0 for a
# column dated t). Lags add up; a difference and a function are dated by their most recent
# input, an interaction by its most recent factor. Columns in the order of `_names`.
_collags(::Union{VarTerm, ConstTerm}) = [0]
function _collags(t::FunTerm)
    ls = [first(_collags(a)) for a in t.args if !(a isa ConstTerm)]
    return [isempty(ls) ? 0 : minimum(ls)]
end
function _collags(t::ProductTerm)
    out = _collags(t.terms[1])
    for c in t.terms[2:end]
        lc = _collags(c)
        out = [min(a, b) for a in out for b in lc]
    end
    return out
end
_collags(t::LagsTerm) = [l + c for l in t.lags for c in _collags(t.term)]
_collags(t::FirstDiffTerm) = _collags(t.term)

# Whether some regressor equals shock column j one period earlier, on the rows where both
# are defined: the shock's first lag, however it is written (`lags(s, 1)`, `lags(s, 0:1)`,
# `lags(firstdiff(s), 1)`, a pre-built column, …).
function _haslag(d::PanelDesign, j::Int)
    xl = _gather(d.X[:, j], _shift(d.pidx, -1))
    for k in axes(d.X, 2)
        k == j && continue
        x = view(d.X, :, k)
        both = 0
        same = true
        for r in eachindex(xl)
            (isfinite(xl[r]) && isfinite(x[r])) || continue
            both += 1
            if !isapprox(xl[r], x[r]; rtol = 1e-12)
                same = false
                break
            end
        end
        same && both > 0 && return true
    end
    return false
end

# Number of groups of a fixed effect whose rows all belong to one unit.
function _singleunitgroups(fcodes::Vector{Int32}, units::AbstractVector{Int32})
    L = Int(maximum(fcodes; init = Int32(0)))
    first_unit = zeros(Int32, L)
    multi = falses(L)
    @inbounds for r in eachindex(fcodes)
        g = fcodes[r]
        if first_unit[g] == 0
            first_unit[g] = units[r]
        elseif first_unit[g] != units[r]
            multi[g] = true
        end
    end
    return count(g -> first_unit[g] != 0 && !multi[g], 1:L)
end

const _CDOLS_HINT = " (or CDOLS(validate = false) to skip the checks)"

# Ugarte-Ruiz (2026): the long difference y_{t+h} − y_{t−1} removes the unit effect, and
# the shock dated t with its lag s_{t−1} leaves the shock innovation. The checks: the
# response is a long difference; each shock is dated right after the base period and its
# first lag is a regressor; the outcome enters the right-hand side only predetermined
# (lagged beyond the base period) and, if the response differences it, only in
# differences (its level carries u_i/(1 − ρ_y)); no fixed effect is nested in the units
# (within estimation would bring the Nickell bias back). That last rule is also checked
# on each horizon's sample, in `_fit_horizon`.
function _validate(est::CDOLS, d::PanelDesign, f)
    est.validate || return nothing
    hint = _CDOLS_HINT
    r = d.response
    (r isa LDiffResponse || r isa CumulResponse) || throw(ArgumentError(
        "CDOLS estimates the long difference y_{t+h} − y_{t−1}: write the response as " *
        "ldiff(y), or cumul(g) for a growth rate g" * hint))
    # the response's own date: ldiff(y) and cumul(g) start at t, ldiff(lags(y, 1)) at t − 1
    ℓr = only(_collags(r.term))
    base = ℓr == 0 ? "t" : "t−$ℓr"
    regs, _, _ = _rewrite_rhs(f.rhs)
    dates = vcat(d.intercept ? [0] : Int[], (_collags(t) for t in regs)...)
    for j in d.shockcols
        nm = d.coefnames[j]
        ℓ = dates[j]
        ℓ == ℓr || throw(ArgumentError(
            "CDOLS: the shock $(repr(nm)) is dated $(ℓ == 0 ? "t" : "t−$ℓ"), but the " *
            "long difference starts at $base: the shock must be dated $base, right " *
            "after the base period. Choose the shock with `shock = …`" *
            (ℓ > ℓr ? ", or lag the response, ldiff(lags(y, $ℓ))" : "") * hint))
        _haslag(d, j) || throw(ArgumentError(
            "CDOLS needs the first lag of the shock $(repr(nm)) among the regressors" *
            (ℓ == 0 ? " (lags($nm, 1))" : "") * ": without it the response is biased " *
            "when the shock is persistent" * hint))
    end
    ys = _vars!(Set{Symbol}(), r.term)
    # outcome variables the response differences (all of ldiff(y); inside cumul, one that
    # enters differenced or at two lags, e.g. cumul(firstdiff(y))) are levels; a growth
    # rate g in cumul(g) is not
    function islevel(v)
        r isa LDiffResponse && return true
        uses = _uses!(Tuple{Int, Bool}[], r.term, Set([v]), 0, false)
        return any(last, uses) || length(unique(first.(uses))) > 1
    end
    for t in regs, v in ys

        uses = _uses!(Tuple{Int, Bool}[], t, Set([v]), 0, false)
        isempty(uses) && continue
        if any(((lag, _),) -> lag <= ℓr, uses)
            cols = [nm for (nm, l) in zip(_names(t), _collags(t)) if l <= ℓr]
            throw(ArgumentError(
                "CDOLS: the outcome enters the right-hand side at or after the response's " *
                "start ($base) in $(join(isempty(cols) ? _names(t) : cols, ", ")); use " *
                "lags of it" * hint))
        end
        islevel(v) && !all(last, uses) &&
            throw(ArgumentError(
                "CDOLS: lagged outcomes enter in differences, lags(firstdiff($v), 1:p); a " *
                "level of $v carries the unit effect ($(join(_names(t), ", ")))" * hint))
    end
    for (k, codes) in enumerate(d.fe)
        rows = findall(!=(0), codes)
        isempty(rows) && continue     # no observations at all: reported by panellp
        fc = codes[rows]
        units = view(d.pidx.unit, rows)
        _isnested(fc, Int(maximum(fc)), units) && throw(ArgumentError(
            "CDOLS does not absorb unit effects (fe($(d.fenames[k])) is nested in the " *
            "units): the long difference removes them, and within estimation would bring " *
            "the Nickell bias back. Time effects are allowed" * hint))
        n1 = _singleunitgroups(fc, units)
        n1 > 0 && @warn "CDOLS: $n1 group(s) of fe($(d.fenames[k])) contain a single " *
              "unit; there the fixed effect is a unit effect, which brings back part of " *
              "the Nickell bias of within estimation"
    end
    return nothing
end

# Cumulative-difference OLS is OLS on the checked design: no unit effects to absorb. A
# fixed effect nested in the units of a horizon's sample (e.g. after `subset` or missing
# values) is a unit effect there.
function _fit_horizon(est::CDOLS, d::PanelDesign, h::Int, y::Vector{Float64},
        mask::BitVector, singletons::Int, absorb::AbsorptionMethod)
    fit = _fit_horizon(WithinFE(), d, h, y, mask, singletons, absorb)
    if est.validate && fit.status === :ok
        for k in axes(fit.dof.nested, 1)
            fit.dof.nested[k, 1] && throw(ArgumentError(
                "CDOLS does not absorb unit effects: in the sample of h = $h, " *
                "fe($(d.fenames[k])) is nested in the units, so within estimation and its " *
                "Nickell bias come back" * _CDOLS_HINT))
        end
    end
    return fit
end

# Almuzara & Sancibrián (2026): the response is a lead of the outcome, whose lags are
# added; the regression has unit effects.
function _validate(::LagAugmented, d::PanelDesign, f)
    d.response isa PlainResponse && throw(ArgumentError(
        "LagAugmented adds lags of the outcome: write the response as leads(y), ldiff(y), " *
        "cumul(g) or anchor(y, z), not a pre-built column"))
    # unit effects, or finer ones nested in the units (e.g. firm × state), whose dummies
    # span the unit dummies
    nested(codes) = (rows = findall(!=(0), codes);
        !isempty(rows) &&
            _isnested(codes[rows], Int(maximum(codes)), view(d.pidx.unit, rows)))
    any(nested, d.fe) || throw(ArgumentError(
        "LagAugmented needs unit fixed effects, fe(<id column>) or a finer fixed effect " *
        "nested in the units (e.g. fe(id) & fe(state)): the Almuzara–Sancibrián " *
        "regression has unit effects"))
    return nothing
end

# Columns of horizon h: the formula's, and the lags up to p(h).
function _activecols(est::_ResolvedLags, d::PanelDesign, h::Int)
    findall(<=(est.lags[h - first(est.horizons) + 1]), d.lagorder)
end

function _fit_horizon(est::_ResolvedLags, d::PanelDesign, h::Int, y::Vector{Float64},
        mask::BitVector, singletons::Int, absorb::AbsorptionMethod)
    return _fit_within(d, h, y, mask, singletons, absorb, _activecols(est, d, h),
        est.lags[h - first(est.horizons) + 1])
end
