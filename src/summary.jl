# Accessors, impulse-response summaries and display.

_fit(r::PanelLocalProjection, h::Integer) = r.fits[_hindex(r.horizons, h)]

function _termindex(names::Vector{String}, term)
    nm = string(term)
    j = findfirst(==(nm), names)
    j === nothing &&
        throw(ArgumentError("no coefficient \"$nm\"; coefficients are $(join(repr.(names), ", "))"))
    return j
end

StatsAPI.coef(r::PanelLocalProjection, h::Integer) = _fit(r, h).coef
StatsAPI.coefnames(r::PanelLocalProjection) = r.coefnames
StatsAPI.nobs(r::PanelLocalProjection, h::Integer) = _fit(r, h).dof.nobs

"Number of observations at each horizon."
nobspath(r::PanelLocalProjection) = [f.dof.nobs for f in r.fits]

"Estimated horizons."
horizons(r::PanelLocalProjection) = r.horizons

"Estimation sample of horizon `h`, as a mask over the data rows."
esample(r::PanelLocalProjection, h::Integer) = _fit(r, h).esample

"""
    coefpath(r; term = first shock) -> Vector{Float64}

The coefficient of `term` across horizons (NaN where it was not estimated).
"""
function coefpath(r::PanelLocalProjection; term = r.shocks[1])
    j = _termindex(r.coefnames, term)
    return [f.coef[j] for f in r.fits]
end

# Standard error from a variance; NaN for a negative one (a two-way variance without PSD
# repair can have negative diagonal elements).
_se(v::Float64) = v >= 0 ? sqrt(v) : NaN

"""
    stderror(cv::PanelLPVcov; term) -> Vector{Float64}

Standard errors of `term` across horizons (NaN where the variance is negative or
undefined).
"""
function StatsAPI.stderror(cv::PanelLPVcov; term)
    j = _termindex(cv.coefnames, term)
    return [_se(V[j, j]) for V in cv.V]
end

function _quantile(df::Float64, level::Float64, dist::Symbol)
    p = 1 - (1 - level) / 2
    (dist === :normal || isinf(df)) && return StatsFuns.norminvcdf(p)
    isnan(df) && return NaN
    return StatsFuns.tdistinvcdf(df, p)
end

"""
    summarize(r, spec_or_vcov; term = first shock, level = 0.95, scale = 1, dist = :auto)

Impulse response of `term`: coefficient, standard error and `level` confidence bounds at
every horizon. `spec_or_vcov` is a variance specification (e.g. `Clustered(:unit)`) or the
result of `vcov`. Critical values come from the specification (t with its degrees of
freedom, or normal); `dist = :normal` forces normal ones (e.g. the 1.96 bands of `xtlp`).
`scale` multiplies coefficients and standard errors.
"""
function summarize(r::PanelLocalProjection, cv::PanelLPVcov; term = r.shocks[1],
        level::Real = 0.95, scale::Real = 1.0, dist::Symbol = :auto)
    _check_option(dist, (:auto, :normal), "dist")
    0 < level < 1 || throw(ArgumentError("level must be in (0, 1)"))
    cv.fingerprint == _fingerprint(r) || throw(ArgumentError(
        "this PanelLPVcov was computed from a different result (other horizons, " *
        "coefficients or sample); compute it with vcov(spec, r) for this one"))
    j = _termindex(r.coefnames, term)
    coef = [scale * f.coef[j] for f in r.fits]
    se = [abs(scale) * _se(V[j, j]) for V in cv.V]
    df = [d[j] for d in cv.df]
    q = [_quantile(d, Float64(level), dist) for d in df]
    return IRFSummary(r.coefnames[j], r.response, _label(cv.spec), Float64(level),
        Float64(scale), collect(r.horizons), coef, se, coef .- q .* se, coef .+ q .* se,
        dist === :normal ? fill(Inf, length(df)) : df, nobspath(r))
end
function summarize(r::PanelLocalProjection, spec::PanelVcov; kw...)
    summarize(r, vcov(spec, r); kw...)
end

Tables.istable(::Type{IRFSummary}) = true
Tables.columnaccess(::Type{IRFSummary}) = true
function Tables.columns(s::IRFSummary)
    (horizon = s.horizon, coef = s.coef, se = s.se,
        lower = s.lower, upper = s.upper, df = s.df, nobs = s.nobs)
end
Tables.columnnames(s::IRFSummary) = (:horizon, :coef, :se, :lower, :upper, :df, :nobs)
Tables.getcolumn(s::IRFSummary, nm::Symbol) = getproperty(Tables.columns(s), nm)
Tables.getcolumn(s::IRFSummary, i::Int) = Tables.columns(s)[i]

_fmt(x::Real) = isnan(x) ? "NaN" : Printf.@sprintf("%.6g", x)
_fmtdf(x::Real) = isinf(x) ? "∞" : isnan(x) ? "NaN" : Printf.@sprintf("%.4g", x)

function Base.show(io::IO, ::MIME"text/plain", s::IRFSummary)
    pct = Printf.@sprintf("%g", 100 * s.level)
    println(io, "Impulse response of $(s.response) to $(s.term)",
        s.scale == 1 ? "" : " (× $(_fmt(s.scale)))")
    println(io, "Variance: $(s.spec); $(pct)% bands")
    hdr = ["h", "coef", "se", "lower", "upper", "df", "nobs"]
    rows = [[string(s.horizon[i]), _fmt(s.coef[i]), _fmt(s.se[i]), _fmt(s.lower[i]),
                _fmt(s.upper[i]), _fmtdf(s.df[i]), string(s.nobs[i])]
            for i in eachindex(s.horizon)]
    w = [maximum(length, [hdr[c]; [row[c] for row in rows]]) for c in eachindex(hdr)]
    println(io, join(lpad.(hdr, w), "  "))
    for (i, row) in enumerate(rows)
        print(io, join(lpad.(row, w), "  "))
        i < length(rows) && println(io)
    end
end

function Base.show(io::IO, ::MIME"text/plain", r::PanelLocalProjection{E}) where {E}
    p = r.panel
    println(io, "Panel local projection ($(nameof(E))), horizons $(r.horizons)")
    println(io, "  response: ", r.response)
    println(io, "  shock:    ", join(r.shocks, ", "))
    println(io, "  regressors: ", join(r.coefnames, ", "))
    isempty(r.fenames) || println(io, "  fixed effects: ", join(r.fenames, ", "))
    println(
        io, "  panel: ", p.nunits, " units (", p.id, "), ", p.ntimes, " periods (", p.time,
        "), ", p.nrows, " rows")
    print(io, "  observations by horizon: ", join(nobspath(r), ", "))
    for note in _fitnotes(r)
        print(io, "\n  note: ", note)
    end
end

function Base.show(io::IO, r::PanelLocalProjection{E}) where {E}
    print(io, "PanelLocalProjection{$(nameof(E))}($(r.response), horizons $(r.horizons))")
end

function Base.show(io::IO, cv::PanelLPVcov)
    print(io, "PanelLPVcov: ", _label(cv.spec), ", horizons ", cv.horizons)
end
