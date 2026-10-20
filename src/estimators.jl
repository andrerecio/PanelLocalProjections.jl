# Per-horizon fits, dispatched on the estimator type.

# A fit without coefficients: no observations (:empty), or a shock that is not identified
# in the horizon's sample (:singular), with `nobs` rows in it.
function _empty_fit(d::PanelDesign, h::Int, mask::BitVector, singletons::Int;
        status::Symbol = :empty, nobs::Int = 0, keep::BitVector = falses(length(d.coefnames)),
        absorbed::BitVector = falses(length(d.coefnames)), iterations::Int = 0,
        converged::Bool = true)
    K = length(d.coefnames)
    G = length(d.groups)
    D = length(d.fe)
    dof = DofInfo(nobs, 0, 0, zeros(Int, D), zeros(Int, D, D), zeros(Int, G),
        Matrix{Bool}(falses(D, G)), singletons)
    stats = ScoreStats(zeros(K, K), [zeros(ng, K) for ng in d.ngroups],
        Tuple{Int, Int}[], Matrix{Float64}[], Symbol[])
    return HorizonFit(h, status, fill(NaN, K), keep, absorbed, fill(NaN, K, K),
        fill(NaN, K), mask, falses(d.ngroups[2]), dof, stats, iterations, converged)
end

# Within transformation, OLS and score statistics on the rows of `mask`.
function _within_fit(d::PanelDesign, y::Vector{Float64}, mask::BitVector,
        absorb::AbsorptionMethod)
    rows = findall(mask)
    Xraw = d.X[rows, :]
    A = hcat(y[rows], Xraw)
    fs = _fe_sample(d.fe, rows)
    it, conv = _absorb_matrix!(A, fs, absorb)
    yt = A[:, 1]
    Xt = A[:, 2:end]
    o = _ols(Xt, yt, Xraw, !isempty(d.fe), d.intercept, d.shockcols)
    return merge(o, (; rows, Xt, yt, fs, xsd = _xsd(Xraw, d.intercept), iterations = it,
        converged = conv))
end

function _singular_fit(d::PanelDesign, h::Int, mask::BitVector, singletons::Int, w;
        status::Symbol = :singular)
    return _empty_fit(d, h, mask, singletons; status, nobs = length(w.rows),
        keep = w.keep, absorbed = w.absorbed, iterations = w.iterations,
        converged = w.converged)
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
function _xsd(Xraw::Matrix{Float64}, intercept::Bool)
    n = size(Xraw, 1)
    sd = [n > 1 ? max(Statistics.std(view(Xraw, :, j)), 1e-3) : 1e-3
          for j in axes(Xraw, 2)]
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
    any(mask) || return _empty_fit(d, h, mask, singletons)
    w = _within_fit(d, y, mask, absorb)
    w.status === :singular && return _singular_fit(d, h, mask, singletons, w)
    U = w.Xt[:, w.keep] .* w.resid
    groups = [view(g, w.rows) for g in d.groups]
    stats, present = _score_stats(U, w.keep, groups, d.ngroups)
    dof = _dofinfo(d, w.rows, w.fs, count(w.keep), present, singletons)
    coef, bread = _expand(w.β, w.B, w.keep)
    return HorizonFit(h, :ok, coef, w.keep, w.absorbed, bread, w.xsd, mask,
        _present(groups[2], d.ngroups[2]), dof, stats, w.iterations, w.converged)
end
