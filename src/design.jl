# Horizon-invariant part of a panel local projection: panel index, numeric columns,
# right-hand side over all data rows, fixed-effect and cluster codes.

function _getcol(cols, name::Symbol)
    name in Tables.columnnames(cols) ||
        throw(ArgumentError("column :$name not found in the data"))
    return Tables.getcolumn(cols, name)
end

_clustername(c::Symbol) = c
_clustername(c::Tuple) = Symbol(join(string.(c), "&"))
_clustervars(c::Symbol) = [c]
_clustervars(c::Tuple) = collect(Symbol, c)

# One cluster specification: a column name, or a tuple of names for their intersection.
_clusterspec(c::Union{Symbol, AbstractString}) = Symbol(c)
function _clusterspec(c::Tuple)
    all(x -> x isa Union{Symbol, AbstractString}, c) && length(c) >= 1 &&
        return Tuple(Symbol.(c))
    throw(ArgumentError("a cluster intersection is a tuple of column names; got $(repr(c))"))
end
function _clusterspec(c)
    throw(ArgumentError(
        "clusters are column names or tuples of column names (an intersection); got $(repr(c))"))
end

# `clusters` is one specification or a collection of them: `:region`, `(:region, :sector)`
# (two dimensions), `((:sector, :quarter),)` (one intersection).
function _clusterlist(clusters)
    clusters isa Union{Symbol, AbstractString} && return Any[Symbol(clusters)]
    clusters isa Union{Tuple, AbstractVector} ||
        throw(ArgumentError("clusters must be a column name or a tuple of cluster " *
                            "specifications; got $(repr(clusters))"))
    return Any[_clusterspec(c) for c in clusters]
end

# Whether two score-dimension names are the same dimension; intersections match whatever
# the order of their columns ("t&grp" is "grp&t").
function _samedim(a::Symbol, b::Symbol)
    a === b && return true
    sa, sb = string(a), string(b)
    (occursin('&', sa) && occursin('&', sb)) || return false
    return sort!(split(sa, '&')) == sort!(split(sb, '&'))
end

function _subsetmask(subset, cols, n)
    subset === nothing && return trues(n)
    v = subset isa Symbol ? _getcol(cols, subset) : subset
    length(v) == n ||
        throw(DimensionMismatch("subset has length $(length(v)); the data have $n rows"))
    eltype(v) <: Union{Missing, Bool} ||
        throw(ArgumentError("subset must be a Bool column or vector"))
    return BitVector([x === true for x in v])
end

function _select_shocks(shock, coefnames, regs)
    if shock === nothing
        isempty(regs) && throw(ArgumentError("the formula has no regressors"))
        length(regs) > 1 && any(t -> t isa ProductTerm, regs) &&
            throw(ArgumentError(
                "the formula has interaction terms, which StatsModels lists after the main " *
                "effects: pass `shock` explicitly (e.g. shock = \"s & X\")"))
        names = _names(first(regs))
    else
        names = shock isa AbstractVector ? string.(shock) : [string(shock)]
        isempty(names) &&
            throw(ArgumentError("`shock` is empty: name at least one coefficient"))
        allunique(names) || throw(ArgumentError(
            "`shock` names a coefficient twice: $(join(repr.(names), ", "))"))
    end
    idx = Int[]
    for nm in names
        j = findfirst(==(nm), coefnames)
        j === nothing && throw(ArgumentError(
            "shock \"$nm\" is not a coefficient; coefficients are $(join(repr.(coefnames), ", "))"))
        push!(idx, j)
    end
    return idx
end

function _build_design(f::StatsModels.FormulaTerm, data; id::Symbol, time::Symbol,
        frequency, subset, clusters, shock)
    cols = Tables.columns(data)
    pidx = PanelIndex(_getcol(cols, id), _getcol(cols, time); frequency = frequency)
    n = nrows(pidx)
    # a column table need not have columns of equal length: check before any lookup
    function column(v::Symbol)
        c = _getcol(cols, v)
        length(c) == n || throw(DimensionMismatch(
            "column :$v has $(length(c)) rows; the id and time columns have $n"))
        return c
    end

    response = _rewrite_lhs(f.lhs)
    regs, fes, ireq, iv = _rewrite_rhs(f.rhs)
    endoterms, instterms = iv === nothing ? (Any[], Any[]) : iv
    _width1(response.term, "the response")
    response isa AnchorResponse && _width1(response.base, "the anchor of anchor(y, z)")

    vars = Set{Symbol}()
    _vars!(vars, response)
    foreach(t -> _vars!(vars, t), regs)
    foreach(t -> _vars!(vars, t), endoterms)
    foreach(t -> _vars!(vars, t), instterms)
    columns = Dict{Symbol, Vector{Float64}}()
    for v in vars
        c = column(v)
        _isnumericcol(c) || throw(ArgumentError(
            "column :$v is not numeric. Categorical regressors are not supported in v0.1: " *
            "absorb them with fe($v) or build dummy columns"))
        columns[v] = _nanfloat(c)
    end
    ctx = _EvalContext(pidx, columns)

    intercept = isempty(fes) ? ireq != :no : false
    ireq === :yes && !isempty(fes) &&
        @warn "the fixed effects absorb the intercept; `1 +` is ignored"
    blocks = Matrix{Float64}[]
    coefnames = String[]
    termof = Tuple{Int, Int}[]               # per column: (term, column within the term)
    if intercept
        push!(blocks, ones(n, 1))
        push!(coefnames, "(Intercept)")
        push!(termof, (0, 1))
    end
    # the endogenous regressors first, as ivreghdfe lists them
    allregs = vcat(endoterms, regs)
    for (k, t) in enumerate(allregs)
        push!(blocks, _cols(t, ctx))
        append!(coefnames, _names(t))
        append!(termof, ((k, c) for c in 1:size(blocks[end], 2)))
    end
    isempty(blocks) && throw(ArgumentError("the formula has no regressors"))
    X = reduce(hcat, blocks)
    allunique(coefnames) ||
        throw(ArgumentError("duplicated regressors: $(coefnames[findall(_nonunique(coefnames))])"))
    nendo = sum(t -> length(_names(t)), endoterms; init = 0)
    endo = collect((1:nendo) .+ intercept)
    Z, instnames = _instruments(instterms, ctx, coefnames, endo)
    shockcols = shock === nothing && !isempty(endo) ? copy(endo) :
                _select_shocks(shock, coefnames, regs)
    any(j -> coefnames[j] == "(Intercept)", shockcols) &&
        throw(ArgumentError("the intercept cannot be the shock"))

    rhsok = BitVector([_rowfinite(X, r) && _rowfinite(Z, r) for r in 1:n])

    # Codes of a fixed-effect or cluster variable: the panel's own unit and period codes
    # for the id and time columns (so that dates group by period, as in leads and lags),
    # the sorted levels of the column otherwise.
    periodcodes = Int32.(pidx.period .- (pidx.pmin - 1))
    function reserved(v::Symbol)
        (v in (:unit, :time) && v !== id && v !== time && v in Tables.columnnames(cols)) ||
            return nothing
        throw(ArgumentError(
            "the column :$v is not the panel's $(v === :unit ? "id" : "time") column, but " *
            ":$v names that panel dimension (e.g. in Clustered(:$v)); rename the column"))
    end
    function varcodes(v::Symbol)
        v === id && return pidx.unit
        v === time && return periodcodes
        reserved(v)
        return _encode(column(v))[1]
    end

    fecodes = Vector{Int32}[]
    fenames = String[]
    for t in fes
        push!(fecodes, _combine_codes([varcodes(v) for v in t.vars])[1])
        push!(fenames, _fename(t))
    end
    allunique(fenames) || throw(ArgumentError("duplicated fixed effects"))

    # Score-sum dimensions: unit and time always; then the FE dimensions and the declared
    # clusters that are not already among them.
    groups = Vector{Int32}[pidx.unit, periodcodes]
    groupnames = Symbol[:unit, :time]
    ngroups = Int[nunits(pidx), nspan(pidx)]
    alias(s::Symbol) = s === id ? :unit : s === time ? :time : s
    stored(nm) = any(g -> _samedim(nm, g), groupnames)
    for (k, t) in enumerate(fes)
        nm = length(t.vars) == 1 ? alias(t.vars[1]) : Symbol(fenames[k])
        stored(nm) && continue
        push!(groups, fecodes[k])
        push!(groupnames, nm)
        push!(ngroups, Int(maximum(fecodes[k]; init = Int32(0))))
    end
    for c in _clusterlist(clusters)
        foreach(reserved, _clustervars(c))
        nm = c isa Symbol ? alias(c) : _clustername(c)
        stored(nm) && continue
        codes, G = _combine_codes([varcodes(v) for v in _clustervars(c)])
        push!(groups, codes)
        push!(groupnames, nm)
        push!(ngroups, G)
    end

    none = (Float64[], Float64[])
    exps = [termof[j][1] == 0 ? none :
            _exposure(allregs[termof[j][1]], termof[j][2], view(X, :, j), ctx)
            for j in axes(X, 2)]
    aggseries, exposures = first.(exps), last.(exps)
    timefe = [k for (k, c) in enumerate(fecodes) if _timenested(c, periodcodes)]

    return PanelDesign(pidx, response, _respname(response), columns, X, coefnames,
        shockcols, intercept, rhsok, fecodes, fenames, groups, groupnames, ngroups,
        _subsetmask(subset, cols, n), aggseries, exposures,
        _owners(coefnames, shockcols, exposures), timefe, zeros(Int, length(coefnames)),
        endo, Z, instnames)
end

# The excluded instruments of an IV block, evaluated on the panel (n × 0 without one), and
# their names. Exogenous regressors instrument themselves: naming one among the excluded
# instruments is an error. An instrument may be an endogenous regressor (that column is
# then exogenous after all: `(x ~ x)` is OLS).
function _instruments(instterms::Vector{Any}, ctx::_EvalContext,
        coefnames::Vector{String}, endo::Vector{Int})
    isempty(instterms) && return zeros(_n(ctx), 0), String[]
    Z = reduce(hcat, [_cols(t, ctx) for t in instterms])
    names = reduce(vcat, [_names(t) for t in instterms])
    allunique(names) ||
        throw(ArgumentError("duplicated instruments: $(names[findall(_nonunique(names))])"))
    for nm in names
        j = findfirst(==(nm), coefnames)
        (j === nothing || j in endo) && continue
        throw(ArgumentError(
            "the instrument $(repr(nm)) is an exogenous regressor; exogenous regressors " *
            "instrument themselves: list only the excluded instruments in (endo ~ inst)"))
    end
    length(names) >= length(endo) || throw(ArgumentError(
        "the IV block has $(length(endo)) endogenous regressor(s) " *
        "($(join(repr.(coefnames[endo]), ", "))) but $(length(names)) excluded " *
        "instrument(s): the model is underidentified"))
    return Z, names
end

# Whether a fixed-effect dimension is nested in the periods (time, or group × time).
function _timenested(codes::Vector{Int32}, periodcodes::Vector{Int32})
    rows = findall(!=(0), codes)
    return _isnested(codes[rows], Int(maximum(codes; init = Int32(0))),
        view(periodcodes, rows))
end

# Whether `x` is an aggregate series: equal across units in every period (where defined).
function _isaggregate(x::AbstractVector{Float64}, pidx::PanelIndex)
    ref = fill(NaN, nspan(pidx))
    @inbounds for r in eachindex(x)
        v = x[r]
        isfinite(v) || continue
        k = pidx.period[r] - pidx.pmin + 1
        a = ref[k]
        if isnan(a)
            ref[k] = v
        elseif v != a && abs(v - a) > 1e-12 * max(abs(v), abs(a))
            return false
        end
    end
    return true
end

# Factor columns of column c of a product with factor widths `w` (first factor slowest,
# as `_rowkron`).
function _kronindex(w::Vector{Int}, c::Int)
    idx = similar(w)
    r = c - 1
    for m in length(w):-1:1
        idx[m] = r % w[m] + 1
        r ÷= w[m]
    end
    return idx
end

# Per calendar period of the span, the value of an aggregate column (equal across units in
# every period); NaN where no row of the period has one.
function _perperiod(x::AbstractVector{Float64}, pidx::PanelIndex)
    a = fill(NaN, nspan(pidx))
    @inbounds for r in eachindex(x)
        v = x[r]
        k = pidx.period[r] - pidx.pmin + 1
        isfinite(v) && isnan(a[k]) && (a[k] = v)
    end
    return a
end

"""
    _exposure(term, c, x, ctx) -> (aggregate series, exposure)

Exposure structure of column `x` (column `c` of `term`): an aggregate series (equal across
units in every period) has exposure 1; an interaction with aggregate factors has the
product of the other factors as exposure, and the product of the aggregate ones as its
aggregate series; a lag of such an interaction lags both. The aggregate series is per
calendar period (`_perperiod`), the exposure per data row. Two empty vectors for any
other column (e.g. a unit-level variable).
"""
function _exposure(t, c::Int, x::AbstractVector{Float64}, ctx::_EvalContext)
    _isaggregate(x, ctx.pidx) && return _perperiod(x, ctx.pidx), ones(length(x))
    return _factorexposure(t, c, ctx, 0)
end

_factorexposure(t, c::Int, ctx::_EvalContext, l::Int) = (Float64[], Float64[])
function _factorexposure(t::LagsTerm, c::Int, ctx::_EvalContext, l::Int)
    k, j = divrem(c - 1, length(_names(t.term)))   # columns: lags outer, inner columns
    return _factorexposure(t.term, j + 1, ctx, l + t.lags[k + 1])
end
function _factorexposure(t::ProductTerm, c::Int, ctx::_EvalContext, l::Int)
    pidx = ctx.pidx
    facs = [_cols(f, ctx) for f in t.terms]
    idx = _kronindex(size.(facs, 2), c)
    shift = l == 0 ? Int32[] : _shift(pidx, -l)
    cols = [l == 0 ? facs[m][:, idx[m]] : _gather(view(facs[m], :, idx[m]), shift)
            for m in eachindex(facs)]
    agg = [_isaggregate(v, pidx) for v in cols]
    any(agg) || return Float64[], Float64[]
    e = ones(nrows(pidx))
    a = ones(nrows(pidx))
    for m in eachindex(cols)
        agg[m] ? (a .*= cols[m]) : (e .*= cols[m])
    end
    return _perperiod(a, pidx), e
end

# Whether two aggregate series are proportional (relative residual 1e-10, on at least two
# periods where both are defined): then, in every period, both columns are their exposure
# times the same number, and the identifying variation of one is its exposure net of the
# other's (Frisch–Waugh–Lovell, period by period).
function _sameagg(a::Vector{Float64}, b::Vector{Float64})
    (isempty(a) || isempty(b)) && return false
    aa = ab = bb = 0.0
    n = 0
    for k in eachindex(a, b)
        (isfinite(a[k]) && isfinite(b[k])) || continue
        aa += a[k]^2
        ab += a[k] * b[k]
        bb += b[k]^2
        n += 1
    end
    (n >= 2 && aa > 0 && bb > 0) || return false
    c = ab / aa
    rr = 0.0
    for k in eachindex(a, b)
        (isfinite(a[k]) && isfinite(b[k])) && (rr += (b[k] - c * a[k])^2)
    end
    return rr <= 1e-20 * bb
end

# Whether two exposures are equal on the rows where both are defined (and there are some).
function _sameexposure(a::Vector{Float64}, b::Vector{Float64})
    n = 0
    @inbounds for r in eachindex(a, b)
        (isfinite(a[r]) && isfinite(b[r])) || continue
        isapprox(a[r], b[r]; rtol = 1e-12) || return false
        n += 1
    end
    return n > 0
end

"""
    _owners(coefnames, shockcols, exposures) -> Vector{Int}

For each coefficient, the column whose synthetic regression its t-LAHR refinement uses:

- a shock, its own;
- a lag of a shock (`Lℓ.shock`), and a regressor with the same exposure as a shock (e.g.
  `s & L1.X` for the shock `s & X`), that shock's;
- any other aggregate series times an exposure (e.g. a control `s2 & X`), its own;
- any other coefficient (outcome lags, unit-level controls), the first shock's.

0 where that column is not an aggregate series times an exposure.
"""
function _owners(coefnames::Vector{String}, shockcols::Vector{Int},
        exposures::Vector{Vector{Float64}})
    has(j) = !isempty(exposures[j])
    owner = zeros(Int, length(coefnames))
    for j in eachindex(coefnames)
        if j in shockcols
            owner[j] = has(j) ? j : 0
            continue
        end
        m = match(r"^L\d+\.(.*)$", coefnames[j])
        k = m === nothing ? nothing :
            findfirst(s -> _wrap(coefnames[s]) == m.captures[1], shockcols)
        if k === nothing && has(j)
            k = findfirst(s -> has(s) && _sameexposure(exposures[j], exposures[s]),
                shockcols)
            k === nothing && (owner[j] = j; continue)
        end
        s = shockcols[something(k, 1)]
        owner[j] = has(s) ? s : 0
    end
    return owner
end

function _nonunique(v::AbstractVector)
    seen = Set{eltype(v)}()
    return [x in seen ? true : (push!(seen, x); false) for x in v]
end
