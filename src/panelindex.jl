# Panel index: unit codes, integer calendar periods and the (unit, period) -> row lookup
# through which every lead, lag and difference is computed. Nothing in the package shifts
# a column by row position.

const _FREQUENCIES = (:year, :halfyear, :quarter, :month)

# Integer period code of a date at the given frequency. Consecutive periods have
# consecutive codes, so leads and lags are integer offsets.
function _period_code(d::Dates.TimeType, frequency::Symbol)
    y = Dates.year(d)
    frequency === :year && return y
    frequency === :halfyear && return 2y + (Dates.month(d) > 6)
    frequency === :quarter && return 4y + Dates.quarterofyear(d) - 1
    frequency === :month && return 12y + Dates.month(d) - 1
    throw(ArgumentError("frequency must be one of $(_FREQUENCIES); got :$frequency"))
end

# First date of the period with the given code (inverse of `_period_code`).
function _period_date(code::Integer, frequency::Symbol)
    frequency === :year && return Dates.Date(code)
    frequency === :halfyear && return Dates.Date(fld(code, 2), 6mod(code, 2) + 1)
    frequency === :quarter && return Dates.Date(fld(code, 4), 3mod(code, 4) + 1)
    frequency === :month && return Dates.Date(fld(code, 12), mod(code, 12) + 1)
    throw(ArgumentError("frequency must be one of $(_FREQUENCIES); got :$frequency"))
end

# Readable label of a period code: 2001, 2001H2, 2001Q3, 2001-07; the code itself for
# integer periods.
function _period_label(code::Integer, frequency)
    frequency === nothing && return string(code)
    frequency === :year && return string(code)
    frequency === :halfyear && return "$(fld(code, 2))H$(mod(code, 2) + 1)"
    frequency === :quarter && return "$(fld(code, 4))Q$(mod(code, 4) + 1)"
    return "$(fld(code, 12))-$(lpad(mod(code, 12) + 1, 2, '0'))"
end

_periods(time::AbstractVector{<:Integer}, frequency) = _periods_integer(time, frequency)
function _periods(time::AbstractVector{<:Dates.TimeType}, frequency)
    frequency === nothing && throw(ArgumentError(
        "the time column holds dates: pass `frequency` (one of $(_FREQUENCIES)) so that " *
        "dates map to consecutive periods"))
    return [_period_code(d, frequency) for d in time]
end
function _periods(time::AbstractVector{<:Real}, frequency)
    all(t -> isfinite(t) && isinteger(t), time) || throw(ArgumentError(
        "the time column must hold integer periods (or dates with `frequency`); found " *
        "non-integral values such as $(first(t for t in time if !(isfinite(t) && isinteger(t))))"))
    return _periods_integer(Int.(time), frequency)
end
function _periods(time::AbstractVector, frequency)
    any(ismissing, time) && throw(ArgumentError(
        "the time column has missing values; drop those rows before calling panellp"))
    v = identity.(time)       # narrow an abstract element type (e.g. Any holding Ints)
    # (each branch dispatches to a method above, so this never recurses)
    (eltype(v) <: Real || eltype(v) <: Dates.TimeType) || throw(ArgumentError(
        "the time column must hold integer periods or Date/DateTime values (with " *
        "`frequency`); got element type $(eltype(v)). Convert strings such as \"2001Q1\" " *
        "to dates, and Dates.Year values to integers, first"))
    return _periods(v, frequency)
end

function _periods_integer(time::AbstractVector{<:Integer}, frequency)
    frequency === nothing || throw(ArgumentError(
        "`frequency` applies only to Date/DateTime time columns; integer periods are " *
        "used as they are (consecutive integers are consecutive periods)"))
    return Vector{Int}(time)
end

const _ROWOF_LIMIT_BYTES = 1 << 30   # 1 GiB

function PanelIndex(id::AbstractVector, time::AbstractVector; frequency = nothing)
    length(id) == length(time) ||
        throw(DimensionMismatch("id and time have different lengths"))
    any(ismissing, id) && throw(ArgumentError(
        "the id column has missing values; drop those rows before calling panellp"))
    frequency === nothing || (frequency = Symbol(frequency))
    frequency === nothing || frequency in _FREQUENCIES ||
        throw(ArgumentError(
            "frequency must be one of $(join(repr.(_FREQUENCIES), ", ")); got $(repr(frequency))"))
    period = _periods(time, frequency)
    labels = sort!(unique(id))
    code = Dict(v => Int32(i) for (i, v) in enumerate(labels))
    unit = Int32[code[v] for v in id]
    n = length(unit)
    n > 0 || throw(ArgumentError("the data have no rows"))
    pmin, pmax = extrema(period)
    nunits = length(labels)
    span = pmax - pmin + 1
    bytes = 4 * nunits * span
    bytes <= _ROWOF_LIMIT_BYTES || throw(ArgumentError(
        "the panel spans $span periods for $nunits units: the period lookup would need " *
        "$(round(bytes / 2^30; digits = 1)) GiB. Is the time column a calendar code " *
        "(e.g. YYYYMM)? Pass dates with `frequency`, or consecutive integer periods."))
    rowof = zeros(Int32, nunits, span)
    dups = String[]
    for r in 1:n
        j = period[r] - pmin + 1
        if rowof[unit[r], j] != 0
            d = "($(repr(id[r])), $(_period_label(period[r], frequency)))"
            length(dups) < 5 && !(d in dups) && push!(dups, d)
        else
            rowof[unit[r], j] = r
        end
    end
    isempty(dups) || throw(ArgumentError(
        "duplicate (id, period) pairs, e.g. $(join(dups, ", ")); each unit can appear " *
        "once per period" *
        (frequency === nothing ? "" : " (dates are grouped into $(frequency) periods)")))
    freq = frequency === nothing ? :none : frequency
    pidx = PanelIndex(unit, period, string.(labels), pmin, pmax, freq, rowof)
    _warn_sparse_span(pidx)
    return pidx
end

nunits(p::PanelIndex) = size(p.rowof, 1)
nspan(p::PanelIndex) = size(p.rowof, 2)
nrows(p::PanelIndex) = length(p.unit)

# Warn when the integer time column looks like a calendar code rather than consecutive
# periods: most of the span is never observed by any unit (YYYYMM codes leave 88 empty
# "months" per year), so every lead or lag across a gap is lost.
function _warn_sparse_span(p::PanelIndex)
    span = nspan(p)
    span <= 1 && return nothing
    seen = falses(span)
    for t in p.period
        seen[t - p.pmin + 1] = true
    end
    unobserved = span - count(seen)
    if unobserved > span ÷ 2
        @warn "$unobserved of the $span periods between $(p.pmin) and $(p.pmax) are " *
              "absent for every unit. If the time column is a calendar code (such as " *
              "YYYYMM), pass a Date column with `frequency` instead, or leads and lags " *
              "across those gaps are undefined."
    end
    return nothing
end

"""
    _shift(p::PanelIndex, k::Integer) -> Vector{Int32}

For each data row (unit i, period t), the data row of (i, t + k), or 0 when that cell is
absent or outside the span. `k > 0` is a lead, `k < 0` a lag.
"""
function _shift(p::PanelIndex, k::Integer)
    out = Vector{Int32}(undef, nrows(p))
    span = nspan(p)
    @inbounds for r in eachindex(out)
        j = p.period[r] - p.pmin + 1 + k
        out[r] = 1 <= j <= span ? p.rowof[p.unit[r], j] : Int32(0)
    end
    return out
end

"""
    _gather(x, idx) -> Vector{Float64}

`x[idx[r]]` for each row, and `NaN` where `idx[r] == 0`.
"""
function _gather(x::AbstractVector{Float64}, idx::AbstractVector{Int32})
    out = Vector{Float64}(undef, length(idx))
    @inbounds for r in eachindex(idx)
        j = idx[r]
        out[r] = j == 0 ? NaN : x[j]
    end
    return out
end

_lead(p::PanelIndex, x::AbstractVector{Float64}, k::Integer) = _gather(x, _shift(p, k))
_lag(p::PanelIndex, x::AbstractVector{Float64}, k::Integer) = _gather(x, _shift(p, -k))

function Base.show(io::IO, p::PanelIndex)
    nobs = nrows(p)
    cells = nunits(p) * nspan(p)
    print(io, "PanelIndex: ", nunits(p), " units × ", nspan(p), " periods (",
        p.frequency === :none ? "$(p.pmin)–$(p.pmax)" :
        "$(_period_date(p.pmin, p.frequency))–$(_period_date(p.pmax, p.frequency)), $(p.frequency)",
        "), ", nobs, " rows",
        nobs == cells ? ", balanced" :
        ", $(cells - nobs) absent cell$(cells - nobs == 1 ? "" : "s")")
end
