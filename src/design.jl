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
    regs, fes, ireq = _rewrite_rhs(f.rhs)
    _width1(response.term, "the response")
    response isa AnchorResponse && _width1(response.base, "the anchor of anchor(y, z)")

    vars = Set{Symbol}()
    _vars!(vars, response)
    foreach(t -> _vars!(vars, t), regs)
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
    if intercept
        push!(blocks, ones(n, 1))
        push!(coefnames, "(Intercept)")
    end
    for t in regs
        push!(blocks, _cols(t, ctx))
        append!(coefnames, _names(t))
    end
    isempty(blocks) && throw(ArgumentError("the formula has no regressors"))
    X = reduce(hcat, blocks)
    allunique(coefnames) ||
        throw(ArgumentError("duplicated regressors: $(coefnames[findall(_nonunique(coefnames))])"))
    shockcols = _select_shocks(shock, coefnames, regs)
    any(j -> coefnames[j] == "(Intercept)", shockcols) &&
        throw(ArgumentError("the intercept cannot be the shock"))

    rhsok = BitVector([_rowfinite(X, r) for r in 1:n])

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

    return PanelDesign(pidx, response, _respname(response), columns, X, coefnames,
        shockcols, intercept, rhsok, fecodes, fenames, groups, groupnames, ngroups,
        _subsetmask(subset, cols, n))
end

function _nonunique(v::AbstractVector)
    seen = Set{eltype(v)}()
    return [x in seen ? true : (push!(seen, x); false) for x in v]
end
