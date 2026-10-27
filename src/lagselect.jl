# Lag augmentation (Almuzara & Sancibrián 2026): the number of lags p(h) of the outcome and
# of every shock added at each horizon, and the information criterion on the synthetic time
# series that selects it. Written from the paper (eq. 15, fn 11; supplement B.1, B.3–B.4,
# B.7); the authors' code (which has no licence) only served as a reference to compare
# outputs with.

function ICLags(; criterion::Symbol = :aic, pmax::Union{Nothing, Integer} = nothing,
        vars::Symbol = :outcome, sample::Symbol = :perp)
    return ICLags(criterion, pmax === nothing ? nothing : Int(pmax), vars, sample)
end

function _label(r::ICLags)
    opts = String[]
    r.criterion === :aic || push!(opts, "criterion = :$(r.criterion)")
    r.pmax === nothing || push!(opts, "pmax = $(r.pmax)")
    r.vars === :outcome || push!(opts, "vars = :$(r.vars)")
    r.sample === :perp || push!(opts, "sample = :$(r.sample)")
    return "ICLags($(join(opts, ", ")))"
end
_label(r::CappedLags) = "CappedLags($(r.pmax === nothing ? "" : r.pmax))"
_label(r::ConstantLags) = "ConstantLags($(r.p))"

# The outcome whose lags are added: y for leads(y) and anchor(y, z), Δy for ldiff(y) (the
# long difference cumulates Δy), g for cumul(g).
_lagvar(r::LeadsResponse) = r.term
_lagvar(r::AnchorResponse) = r.term
_lagvar(r::CumulResponse) = r.term
_lagvar(r::LDiffResponse) = FirstDiffTerm(r.term)

# ⌈(T − H)^{1/3}⌉, T the distinct periods in the data and H the largest horizon: the
# default p_max of the Almuzara–Sancibrián code.
function _code_pmax(d::PanelDesign, horizons::UnitRange{Int})
    T = length(unique(d.pidx.period))
    return ceil(Int, max(T - last(horizons), 0)^(1 / 3))
end

# Whether column c equals x on every row where x is defined (and x is defined somewhere).
function _samecol(x::AbstractVector{Float64}, c::AbstractVector{Float64})
    n = 0
    @inbounds for r in eachindex(x, c)
        isfinite(x[r]) || continue
        (isfinite(c[r]) && isapprox(x[r], c[r]; rtol = 1e-12)) || return false
        n += 1
    end
    return n > 0
end

# The aggregate series a, lagged l periods.
function _lagseries(a::Vector{Float64}, l::Int)
    out = fill(NaN, length(a))
    out[(l + 1):end] = a[1:(end - l)]
    return out
end

"""
    _with_lags(d, ctx, pmax) -> PanelDesign

The design with lags 1:pmax of the outcome (`_lagvar`) and of every shock column
appended, named `L1.y`, `L1.(s & X)`, … (`lagorder` = the lag). A lag already among the
regressors, whatever its name (e.g. `s & L1.X` for the lag of the shock `s & L0.X`), is not
added again. The lag of an aggregate series times an exposure lags both.

Panel LP-IV (Almuzara & Sancibrián 2026, §3.4): the lags of an endogenous shock are
endogenous, and for each lag ℓ of an endogenous shock that is added, lag ℓ of every excluded
instrument is appended to the instruments (`zlagorder`), unless equal to one of them; the
lags of the outcome are exogenous. So the fit at horizon h depends on p(h) only.
"""
function _with_lags(d::PanelDesign, ctx::_EvalContext, pmax::Int)
    pmax == 0 && return d
    lv = _lagvar(d.response)
    vars = [(only(_names(lv)), _col1(lv, ctx), Float64[], Float64[], false)]
    for j in d.shockcols
        push!(vars, (
            d.coefnames[j], d.X[:, j], d.exposures[j], d.aggseries[j], j in d.endo))
    end
    cols = Vector{Float64}[]
    names = String[]
    order = Int[]
    lagexp = Vector{Float64}[]
    lagagg = Vector{Float64}[]
    isendo = Bool[]
    for (nm, v, e, a, en) in vars, l in 1:pmax

        lname = string("L", l, ".", _wrap(nm))
        (lname in d.coefnames || lname in names) && continue
        shift = _shift(d.pidx, -l)
        x = _gather(v, shift)
        (any(j -> _samecol(x, view(d.X, :, j)), axes(d.X, 2)) ||
         any(c -> _samecol(x, c), cols)) && continue
        push!(cols, x)
        push!(names, lname)
        push!(order, l)
        push!(lagexp, isempty(e) ? Float64[] : _gather(e, shift))
        push!(lagagg, isempty(a) ? Float64[] : _lagseries(a, l))
        push!(isendo, en)
    end
    Z, instnames, zlag = d.Z, d.instnames, d.zlagorder
    zagg, zexp = d.zaggseries, d.zexposures
    if any(isendo)
        zc = Vector{Float64}[]
        endolags = Set(order[isendo])
        for k in eachindex(d.instnames), l in 1:pmax

            (d.zlagorder[k] == 0 && l in endolags) || continue
            shift = _shift(d.pidx, -l)
            x = _gather(view(d.Z, :, k), shift)
            (any(c -> _samecol(x, view(d.Z, :, c)), axes(d.Z, 2)) ||
             any(c -> _samecol(x, c), zc)) && continue
            push!(zc, x)
            instnames = vcat(instnames, string("L", l, ".", _wrap(d.instnames[k])))
            zlag = vcat(zlag, l)
            e, a = d.zexposures[k], d.zaggseries[k]
            zexp = vcat(zexp, [isempty(e) ? Float64[] : _gather(e, shift)])
            zagg = vcat(zagg, [isempty(a) ? Float64[] : _lagseries(a, l)])
        end
        isempty(zc) || (Z = hcat(d.Z, reduce(hcat, zc)))
    end
    isempty(cols) && return d
    K = length(d.coefnames)
    coefnames = vcat(d.coefnames, names)
    exposures = vcat(d.exposures, lagexp)
    endo = vcat(d.endo, K .+ findall(isendo))
    return PanelDesign(d.pidx, d.response, d.responsename, d.columns,
        hcat(d.X, reduce(hcat, cols)), coefnames, d.shockcols, d.intercept, d.rhsok, d.fe,
        d.fenames, d.groups, d.groupnames, d.ngroups, d.subset, vcat(d.aggseries, lagagg),
        exposures, _owners(coefnames, d.shockcols, exposures), d.timefe,
        vcat(d.lagorder, order), endo, Z, instnames, zlag, zagg, zexp)
end

"""
    _reducedform(d) -> PanelDesign

The reduced form of a panel LP-IV design: the excluded instruments in place of the
endogenous columns, right after the intercept (if any). The endogenous shocks are replaced,
at the place of the first of them among the shocks, by the first instrument with the same
exposure as that shock (`_sameexposure`: `s & mp` for `s & ffr`), which therefore weights
the `ICLags` criterion when it leads; an error when there is none. The other instruments
are plain regressors (one absorbed by the fixed effects is omitted, as in the IV fit). The
design itself when it has no IV block.
"""
function _reducedform(d::PanelDesign)
    isempty(d.endo) && return d
    lead = d.intercept ? [1] : Int[]
    exo = [j for j in axes(d.X, 2) if !(j in d.endo) && !(j in lead)]
    L = size(d.Z, 2)
    X = hcat(d.X[:, lead], d.Z, d.X[:, exo])
    coefnames = vcat(d.coefnames[lead], d.instnames, d.coefnames[exo])
    off = length(lead)
    pos(j) = off + L + findfirst(==(j), exo)
    shockcols = Int[]
    for j in d.shockcols
        if !(j in d.endo)
            push!(shockcols, pos(j))
        elseif !any(k -> k <= off + L, shockcols)
            e1 = d.exposures[j]
            k1 = isempty(e1) ? nothing :
                 findfirst(
                k -> !isempty(d.zexposures[k]) && _sameexposure(d.zexposures[k], e1),
                1:L)
            k1 === nothing && throw(ArgumentError(
                "ICLags on panel LP-IV weights the reduced form with the exposure of the " *
                "shock $(repr(d.coefnames[j])), and no excluded instrument has that exposure " *
                "(an aggregate series times the same unit-level variable): use ConstantLags " *
                "or CappedLags"))
            push!(shockcols, off + k1)
        end
    end
    exposures = vcat(d.exposures[lead], d.zexposures, d.exposures[exo])
    return PanelDesign(d.pidx, d.response, d.responsename, d.columns, X, coefnames,
        shockcols, d.intercept, d.rhsok, d.fe, d.fenames, d.groups, d.groupnames,
        d.ngroups, d.subset, vcat(d.aggseries[lead], d.zaggseries, d.aggseries[exo]),
        exposures, _owners(coefnames, shockcols, exposures), d.timefe,
        vcat(d.lagorder[lead], d.zlagorder, d.lagorder[exo]), Int[],
        zeros(size(d.X, 1), 0), String[], Int[], Vector{Float64}[], Vector{Float64}[])
end

# `mask` restricted to the rows where the lags of horizon h are defined.
_lagmask!(mask::BitVector, ::PanelLPEstimator, d::PanelDesign, h::Int) = mask
function _lagmask!(mask::BitVector, est::_ResolvedLags, d::PanelDesign, h::Int)
    p = est.lags[h - first(est.horizons) + 1]
    for j in eachindex(d.lagorder)
        0 < d.lagorder[j] <= p || continue
        mask .&= isfinite.(view(d.X, :, j))
    end
    for k in eachindex(d.zlagorder)
        0 < d.zlagorder[k] <= p || continue
        mask .&= isfinite.(view(d.Z, :, k))
    end
    return mask
end

"""
    _augment(estimator, d, ctx, base, horizons, singletons, absorb)
        -> (design, estimator, LagSelection or nothing)

For `LagAugmented`: resolve p(h), and append the lag columns to the design. Other
estimators pass through.
"""
function _augment(
        est::PanelLPEstimator, d::PanelDesign, ctx, base, horizons, singletons, absorb)
    (d, est, nothing)
end
function _augment(est::LagAugmented, d::PanelDesign, ctx::_EvalContext, base::BitVector,
        horizons::UnitRange{Int}, singletons::Symbol, absorb::AbsorptionMethod)
    # panel LP-IV: the information criterion works on the reduced form (the instruments in
    # place of the endogenous regressors); the other rules need no regression
    dl = est.rule isa ICLags ? _reducedform(d) : d
    sel, lags = _lagrule(est.rule, dl, ctx, base, horizons, singletons, absorb)
    return _with_lags(d, ctx, maximum(lags)), _ResolvedLags(lags, horizons), sel
end

function _lagrule(r::ConstantLags, d, ctx, base, horizons, singletons, absorb)
    return LagSelection(_label(r), r.p, r.p, Float64[], Float64[], Int[]),
    fill(r.p, length(horizons))
end

function _lagrule(r::CappedLags, d, ctx, base, horizons, singletons, absorb)
    pmax = something(r.pmax, _code_pmax(d, horizons))
    return LagSelection(_label(r), pmax, pmax, Float64[], Float64[], Int[]),
    [min(h, pmax) for h in horizons]
end

"""
    _lagrule(r::ICLags, …) -> (LagSelection, p(h) per horizon)

Almuzara & Sancibrián (2026, fn 11; supplement B.3–B.4): for p = 0, …, pmax, the h = 0
regression of the outcome on the formula's right-hand side and p of its lags (and of the
shocks, `vars = :outcome_and_shock`), with the formula's fixed effects, on the rows where
they are defined (`sample = :perp`) or on the pmax-lag rows (`:common`); its residuals
averaged per period with the first shock's synthetic weights, ξ̄_t, and
σ̂²(p) = Σ_t ω_t ξ̄_t² / Σ_t ω_t over the periods with positive weight ω_t (fn 24: the
synthetic regression is weighted by ω_t; the mean of ξ̄_t² on balanced panels);
IC(p) = log σ̂²(p) + k(p)·g, with k(p) the lags candidate p adds to the regression (lags
already among the regressors, or collinear with them, add none) and g = 2/T (AIC) or
log(T)/T (BIC), T the periods of the p = 0 regression. A candidate that adds no lag to the
previous one is its regression, with its criterion; p̂ is the last candidate with the
minimising regression (so lags already in the formula count towards it), and
p(h) = max{1, min(h, p̂)} (eq. 15).
"""
function _lagrule(r::ICLags, d::PanelDesign, ctx::_EvalContext, base::BitVector,
        horizons::UnitRange{Int}, singletons::Symbol, absorb::AbsorptionMethod)
    j1 = d.shockcols[1]
    isempty(d.exposures[j1]) && throw(ArgumentError(
        "ICLags averages the residuals with the exposure of the first shock, " *
        "$(repr(d.coefnames[j1])), which is not an aggregate series times a " *
        "unit-level exposure: use ConstantLags or CappedLags"))
    pmax = something(r.pmax, max(1, _code_pmax(d, horizons)))
    da = _with_lags(d, ctx, pmax)
    lv = _lagvar(d.response)
    v = _col1(lv, ctx)
    vname = only(_names(lv))
    # the added lags a candidate counts: of the outcome (and of the shocks)
    counted(j) = da.lagorder[j] > 0 &&
                 (r.vars === :outcome_and_shock ||
                  da.coefnames[j] == string("L", da.lagorder[j], ".", _wrap(vname)))
    # candidate p: the formula's columns, and the counted lags up to p
    candidate(p) = [j
                    for j in eachindex(da.lagorder)
                    if da.lagorder[j] == 0 || (da.lagorder[j] <= p && counted(j))]
    rowsok(p) = reduce(
        .&, (isfinite.(view(da.X, :, j)) for j in candidate(p)
        if da.lagorder[j] > 0);
        init = base .& isfinite.(v))
    common = r.sample === :common ? rowsok(pmax) : nothing
    sigma2 = fill(NaN, pmax + 1)
    periods = zeros(Int, pmax + 1)
    nlags = zeros(Int, pmax + 1)
    model = collect(0:pmax)          # the first candidate with the same regression
    num, den = zeros(nspan(da.pidx)), zeros(nspan(da.pidx))
    prev = nothing
    for p in 0:pmax
        mask = common === nothing ? rowsok(p) : copy(common)
        singletons === :drop && _drop_singletons!(mask, d.fe)
        if !any(mask)
            p == 0 || continue
            why = r.sample === :common ?
                  "ICLags(sample = :common): no row has the formula's variables and " *
                  "$pmax lags of $vname defined; lower pmax" :
                  "ICLags: the h = 0 regression of $vname on the formula without lags " *
                  "has no observations"
            throw(ArgumentError(why))
        end
        w = _within_fit(da, v, mask, absorb, candidate(p))
        if w.status !== :ok
            p == 0 && _icsingular(da, w, vname)
            continue
        end
        nlags[p + 1] = count(j -> w.keep[j] && counted(j), eachindex(w.keep))
        if prev !== nothing && prev.mask == mask && prev.keep == w.keep
            # no lag of candidate p entered the regression: candidate prev.p's model
            sigma2[p + 1], periods[p + 1] = sigma2[prev.p + 1], periods[prev.p + 1]
            model[p + 1] = model[prev.p + 1]
            continue
        end
        prev = (; p, mask, w.keep)
        wt, _ = _shockweights(da, j1, w.rows, w.keep)
        fill!(num, 0.0)
        fill!(den, 0.0)
        for (i, row) in enumerate(w.rows)
            k = da.pidx.period[row] - da.pidx.pmin + 1
            num[k] += wt[i] * w.resid[i]
            den[k] += wt[i]^2
        end
        ts = findall(>(0), den)
        isempty(ts) && continue
        periods[p + 1] = length(ts)
        sigma2[p + 1] = sum(num[t]^2 / den[t] for t in ts) / sum(den[t] for t in ts)
    end
    T = periods[1]
    T > 0 || throw(ArgumentError(
        "ICLags: in the h = 0 regression of $vname on the formula without lags, no period " *
        "has a positive synthetic weight (the exposure of $(repr(da.coefnames[j1])) does " *
        "not vary within any period)"))
    g = r.criterion === :aic ? 2 / T : log(T) / T
    ic = [isnan(s) ? Inf : log(s) + k * g for (k, s) in zip(nlags, sigma2)]
    # the last candidate with the selected regression: lags already in the formula count
    phat = findlast(==(model[argmin(ic)]), model) - 1
    sel = LagSelection(_label(r), pmax, phat, ic, sigma2, periods)
    return sel, [max(1, min(h, phat)) for h in horizons]
end

# The error of an h = 0 criterion regression that does not identify the shocks, worded as
# `_check_identified`'s.
function _icsingular(d::PanelDesign, w, vname::String)
    absorbed = [d.coefnames[j] for j in d.shockcols if w.absorbed[j]]
    why = isempty(absorbed) ?
          "the shock $(join(repr.(d.coefnames[d.shockcols]), ", ")) is collinear with the " *
          "other regressors" :
          "the shock $(join(repr.(absorbed), ", ")) is absorbed by the fixed effects" *
          (isempty(d.timefe) ? "" :
           " (an aggregate shock with time effects? Interact it with a unit-level " *
           "exposure, or drop the time effects)")
    throw(ArgumentError(
        "ICLags: in the h = 0 regression of $vname on the formula without lags, $why"))
end

"""
    lagorder(r) -> Vector{Int}

The number of lags of the outcome and of each shock added at each horizon by
[`LagAugmented`](@ref) (zeros for the other estimators). The rule and, for
[`ICLags`](@ref), the criterion of each candidate are in `r.lags`, a [`LagSelection`](@ref).
"""
lagorder(r::PanelLocalProjection) = [f.nlags for f in r.fits]
