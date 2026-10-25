"""
    panellp(f, data; id, time, horizon, estimator = WithinFE(), shock = nothing,
            frequency = nothing, subset = nothing, sample = :horizon, singletons = :drop,
            absorb = AlternatingProjections(), clusters = ())

Estimate a panel local projection: for each horizon `h`, regress the response defined by
the left-hand side of `f` on its right-hand side, absorbing the `fe()` terms.

# Formula
- Left-hand side: `leads(y)` (y_{i,t+h}), `ldiff(y)` (y_{i,t+h} − y_{i,t−1}),
  `anchor(y, z)` (y_{i,t+h} − z_{i,t}), `cumul(y)` (Σ_{k≤h} y_{i,t+k}), or a plain column
  (a pre-built outcome) with `horizon = 0`.
- Right-hand side: columns and functions of columns, interactions `s & X`, `lags(x, p)`,
  `lags(x, a:b)`, `firstdiff(x)`, and fixed effects `fe(v)`, `fe(a) & fe(b)`. An intercept
  is included when there are no fixed effects (unless `0 +`).
- Panel LP-IV: one block `(endo ~ inst)`, e.g. `(s & ffr ~ s & mp)`, makes the fit two-stage
  least squares (Almuzara & Sancibrián 2026, §3.4) on each horizon's fixed-effects
  design: the endogenous regressors on the left (listed first among the coefficients,
  and the default `shock`), the excluded instruments on the right; the other regressors
  are exogenous and instrument themselves. Every variance applies to the IV scores
  X̂ ξ̂ (ξ̂ the structural residual); `TLAHR(refine = false)` is the paper's LP-IV variance.
  An instrument absorbed by the fixed effects or collinear with the exogenous regressors
  (or with earlier instruments) is left out of the first stage, with a note; a horizon with
  fewer instruments than endogenous regressors left is `:underidentified`. Only with
  `WithinFE`.

# Keywords
- `id`, `time`: the unit and period columns. Integer periods are consecutive integers;
  `Date`/`DateTime` periods need `frequency` (`:year`, `:halfyear`, `:quarter`, `:month`).
  `fe(time)` and clusters on `time` group by period, not by date.
- `horizon`: `H` (horizons `0:H`) or a range of consecutive horizons such as `1:6`.
- `estimator`: [`WithinFE`](@ref) (fixed effects; pooled OLS without them),
  [`SplitPanelJackknife`](@ref), [`CDOLS`](@ref) or [`LagAugmented`](@ref).
- `shock`: coefficient name(s) of the shock; by default the columns of the first
  regressor (the endogenous ones with an IV block). Required when the formula has
  interaction terms and other regressors.
- `subset`: a Bool column or vector restricting the estimation sample *after* leads and
  lags are computed (they may read outside it).
- `sample`: `:horizon` (each horizon its own complete cases) or `:common` (the
  intersection over horizons).
- `singletons`: `:drop` (iteratively, as `reghdfe` and `xtlp`) or `:keep`.
- `absorb`: [`AlternatingProjections`](@ref) or [`OneSweep`](@ref).
- `clusters`: extra columns (or tuples of columns, for intersections) whose score sums are
  stored for clustered variances, e.g. `clusters = (:region, (:sector, :quarter))`; unit,
  time and the fixed-effect dimensions always are.

Leads, lags and differences are taken on the (unit, calendar period) grid: an absent
period makes them undefined, and rows where the response or any regressor is undefined
(missing, `NaN`, or a missing input of a function term) leave that horizon's sample.

Horizons without observations, horizons where the shock is absorbed by the fixed effects
or collinear in the horizon's sample (status `:singular`), omitted regressors and
non-convergence are reported in one warning each per call, and their coefficients are
NaN; it is an error when no horizon identifies the shock.
"""
function panellp(f::StatsModels.FormulaTerm, data; id::Symbol, time::Symbol, horizon,
        estimator::PanelLPEstimator = WithinFE(), shock = nothing, frequency = nothing,
        subset = nothing, sample::Symbol = :horizon, singletons::Symbol = :drop,
        absorb::AbsorptionMethod = AlternatingProjections(), clusters = ())
    _check_option(sample, (:horizon, :common), "sample")
    _check_option(singletons, (:drop, :keep), "singletons")
    (id === :time || time === :unit || id === time) && throw(ArgumentError(
        "`:unit` and `:time` name the panel dimensions (e.g. in Clustered(:unit, :time)): " *
        "the id column cannot be called `time`, nor the time column `unit`"))
    horizons = _horizons(horizon)
    d = _build_design(f, data; id, time, frequency, subset, clusters, shock)
    isempty(d.endo) || _check_iv(estimator)
    _validate(estimator, d, f)
    d.response isa PlainResponse && horizons != 0:0 &&
        throw(ArgumentError(
            "a plain response is a pre-built outcome and is estimated at horizon 0 only; use " *
            "leads(y), ldiff(y), anchor(y, z) or cumul(y) to estimate several horizons"))
    ctx = _EvalContext(d.pidx, d.columns)
    base = d.rhsok .& d.subset
    for codes in Iterators.flatten((d.fe, d.groups))
        base .&= codes .!= 0
    end
    # LagAugmented: p(h) and the lag columns (the other estimators pass through)
    d, est, lags = _augment(estimator, d, ctx, base, horizons, singletons, absorb)
    # the responses are built again when each horizon is fitted, not kept for all
    masks = [_lagmask!(base .& isfinite.(_lhs_at(d.response, ctx, h)), est, d, h)
             for h in horizons]
    if sample === :common
        common = _common_sample(masks, horizons)
        masks = [copy(common) for _ in masks]
    end
    fits = _estimate(est, d, ctx, horizons, masks, singletons, absorb)
    _check_identified(fits, d, time)
    p = d.pidx
    panel = PanelSummary(
        id, time, p.frequency, nunits(p), nspan(p), length(unique(p.period)),
        nrows(p), p.pmin, p.pmax)
    r = PanelLocalProjection(estimator, absorb, f, d.responsename, d.coefnames,
        d.intercept, d.coefnames[d.shockcols], horizons, fits, d.groupnames, d.ngroups,
        d.fenames, panel, sample, singletons, lags, d.owner, d.lagorder,
        d.coefnames[d.endo], d.instnames)
    foreach(note -> @warn(note), _fitnotes(r))
    return r
end

function _horizons(horizon)
    if horizon isa Integer
        horizon >= 0 || throw(ArgumentError("horizon must be non-negative; got $horizon"))
        return 0:Int(horizon)
    elseif horizon isa AbstractVector{<:Integer} && !isempty(horizon)
        hs = collect(Int, horizon)
        a, b = first(hs), last(hs)
        (a >= 0 && hs == a:b) && return a:b
    end
    throw(ArgumentError(
        "horizon must be an integer H ≥ 0 (horizons 0:H) or a range of consecutive " *
        "horizons a:b with a ≥ 0; got $(repr(horizon))"))
end

# Intersection of the per-horizon samples; an error naming the horizon that empties it.
function _common_sample(masks::Vector{BitVector}, horizons::UnitRange{Int})
    common = copy(masks[1])
    for (k, m) in enumerate(masks)
        common .&= m
        any(common) && continue
        what = any(m) ?
               "the intersection of the samples of horizons " *
               "$(first(horizons))–$(horizons[k]) is empty" :
               "horizon $(horizons[k]) has no observations"
        throw(ArgumentError(
            "sample = :common: $what, so no row has the response and every regressor " *
            "defined at all horizons. Use fewer horizons or sample = :horizon."))
    end
    return common
end

# Function barrier: concrete types from here on.
function _estimate(est::PanelLPEstimator, d::PanelDesign, ctx::_EvalContext,
        horizons::UnitRange{Int}, masks::Vector{BitVector}, singletons::Symbol,
        absorb::AbsorptionMethod)
    fits = Vector{HorizonFit}(undef, length(horizons))
    for (k, h) in enumerate(horizons)
        mask = copy(masks[k])
        dropped = singletons === :drop ? _drop_singletons!(mask, d.fe) : 0
        fits[k] = _fit_horizon(est, d, h, _lhs_at(d.response, ctx, h), mask, dropped,
            absorb)
    end
    return fits
end

# An error when no horizon has observations, or when the shock is identified at none.
function _check_identified(fits::Vector{HorizonFit}, d::PanelDesign, time::Symbol)
    all(f -> f.status === :empty, fits) && throw(ArgumentError(
        "no horizon has observations: at every row the response or a regressor is " *
        "undefined, or only singletons are left (check `horizon`, `subset`, the lags and " *
        "missing values)"))
    any(f -> f.status === :ok, fits) && return nothing
    any(f -> f.status === :spj_half_singular, fits) && throw(ArgumentError(
        "no horizon can be estimated: at every horizon the shock is not identified or a " *
        "split-panel half is rank-deficient ($(_halfcause(fits, d.coefnames)))"))
    any(f -> f.status === :singular, fits) || throw(ArgumentError(
        "the model is underidentified at every horizon: once the fixed effects and the " *
        "exogenous regressors are partialled out, fewer excluded instruments than " *
        "endogenous regressors are left, or the instruments do not predict every " *
        "endogenous regressor (instruments $(join(repr.(d.instnames), ", ")))"))
    f = fits[findfirst(f -> f.status === :singular, fits)]
    absorbed = [d.coefnames[j] for j in d.shockcols if f.absorbed[j]]
    if isempty(absorbed)
        throw(ArgumentError(
            "the shock $(join(repr.(d.coefnames[d.shockcols]), ", ")) is collinear with the " *
            "other regressors at every horizon"))
    end
    timefe = any(nm -> string(time) in split(nm, '&'), d.fenames)
    throw(ArgumentError(
        "the shock $(join(repr.(absorbed), ", ")) is absorbed by the " *
        "$(isempty(d.fe) ? "constant" : "fixed effects") at every horizon" *
        (timefe ?
         " (an aggregate shock with time effects? Interact it with a unit-level " *
         "exposure, or drop the time effects)" : "")))
end

# Why split-panel halves are rank-deficient: the regressors the fixed effects absorb inside
# a half (marked in the fits' `absorbed`), or else too few rows.
function _halfcause(fits::Vector{HorizonFit}, coefnames::Vector{String})
    hf = [f for f in fits if f.status === :spj_half_singular]
    names = unique(coefnames[j] for f in hf for j in findall(f.absorbed .& f.keep))
    isempty(names) &&
        return "fewer rows than regressors, or collinear regressors, inside a half: too " *
               "few periods per unit?"
    # every regressor absorbed: units with a single row in a half
    allabs = any(f -> any(f.keep) && all(f.absorbed[f.keep]), hf)
    return "$(join(names, ", ")) absorbed by the fixed effects inside a half" *
           (allabs ? ": too few periods per unit?" : "")
end

# "0–3, 5, 7–9"
function _hlist(hs::Vector{Int})
    parts = String[]
    k = 1
    while k <= length(hs)
        j = k
        while j < length(hs) && hs[j + 1] == hs[j] + 1
            j += 1
        end
        push!(parts, j == k ? string(hs[k]) :
                     j == k + 1 ? "$(hs[k]), $(hs[j])" :
                     "$(hs[k])–$(hs[j])")
        k = j + 1
    end
    return join(parts, ", ")
end

# What went wrong at some horizons: reported by `panellp` (one warning each) and `show`.
function _fitnotes(r::PanelLocalProjection)
    notes = String[]
    hs(pred) = [f.h for f in r.fits if pred(f)]
    estimated = hs(f -> f.status === :ok)
    h = hs(f -> f.status === :empty)
    isempty(h) || push!(notes,
        "no observations at h = $(_hlist(h)) (the response or a regressor is undefined " *
        "at every row, or only singletons are left); coefficients are NaN")
    h = hs(f -> f.status === :singular)
    isempty(h) || push!(notes,
        "the shock is absorbed by the fixed effects or collinear in the sample of " *
        "h = $(_hlist(h)); coefficients are NaN")
    h = hs(f -> f.status === :underidentified)
    isempty(h) || push!(notes,
        "underidentified at h = $(_hlist(h)) (once the fixed effects and the exogenous " *
        "regressors are partialled out, fewer excluded instruments than endogenous " *
        "regressors are left, or the instruments do not predict every endogenous " *
        "regressor); coefficients are NaN")
    h = hs(f -> f.status === :spj_half_singular)
    isempty(h) || push!(notes,
        "a split-panel half is rank-deficient at h = $(_hlist(h)) " *
        "($(_halfcause(r.fits, r.coefnames))); coefficients are NaN")
    absorbed = (f, j) -> f.absorbed[j]
    # a lag LagAugmented added that the other regressors span is redundant, not omitted
    collinear = (f, j) -> f.active[j] && !f.keep[j] && !f.absorbed[j] && r.addedlag[j] == 0
    for (label, pred) in (("absorbed by the fixed effects or the constant", absorbed),
        ("collinear with the other regressors", collinear))
        items = String[]
        for (j, nm) in enumerate(r.coefnames)
            h = hs(f -> f.status === :ok && pred(f, j))
            isempty(h) && continue
            push!(items, h == estimated ? nm : "$nm (h = $(_hlist(h)))")
        end
        isempty(items) || push!(notes, "omitted, $label: $(join(items, "; "))")
    end
    items = String[]
    for (l, nm) in enumerate(r.instruments)
        h = hs(f -> f.status === :ok && !f.instkeep[l])
        isempty(h) && continue
        push!(items, h == estimated ? nm : "$nm (h = $(_hlist(h)))")
    end
    isempty(items) || push!(notes,
        "instruments left out of the first stage (absorbed by the fixed effects, or " *
        "collinear with the exogenous regressors or the other instruments): " *
        join(items, "; "))
    h = hs(f -> !f.converged)
    isempty(h) || push!(notes,
        "fixed-effect absorption did not converge at h = $(_hlist(h)); raise `maxiter` or " *
        "loosen `tol` in AlternatingProjections")
    return notes
end
