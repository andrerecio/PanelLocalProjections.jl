# Formula interface. `@formula` (StatsModels) parses the expression; `_rewrite_lhs` and
# `_rewrite_rhs` turn its terms into the internal tree below, and `_cols` evaluates that
# tree on the panel. Every shift goes through the PanelIndex (`_shift`/`_gather`), so leads
# and lags respect gaps and never cross units.
#
# The left-hand-side helpers (`leads`, `ldiff`, `anchor`, `cumul`) and `firstdiff` follow
# LocalProjections.jl (gragusa; andrerecio fork, src/terms.jl), where they act on a single
# time series; here they act on the (unit, period) grid.

# ---------------------------------------------------------------------------------------
# User-facing helpers. Inside `@formula` they are never called: the macro records them as
# StatsModels.FunctionTerm{typeof(helper)} and `_rewrite_*` interprets them. Called
# directly they build the same terms, which are StatsModels terms, so formulas can be
# built programmatically:
#     leads(term(:y)) ~ term(:x) + lags(term(:x), p) + fe(term(:id))
# ---------------------------------------------------------------------------------------

_term(s::Symbol) = VarTerm(s)
_term(t::StatsModels.AbstractTerm) = _inner(t)
_term(t) = t

"""
    leads(y)

Response `y_{i,t+h}` at horizon `h` (left-hand side of a `panellp` formula).
"""
leads(y) = LeadsResponse(_term(y))

"""
    ldiff(y)

Long-difference response `y_{i,t+h} − y_{i,t−1}` at horizon `h` (Stata `xtlp`
`ytransf(cmltdiff)`, `locproj` `tr(cmlt)`).
"""
ldiff(y) = LDiffResponse(_term(y))

"""
    anchor(y, z)

Anchored response `y_{i,t+h} − z_{i,t}`; `anchor(y, y)` is the h-difference
`y_{i,t+h} − y_{i,t}`.
"""
anchor(y, z) = AnchorResponse(_term(y), _term(z))

"""
    cumul(y)

Cumulated response `Σ_{k=0}^{h} y_{i,t+k}` (Stata `xtlp` `ytransf(cmltsum)`); with a growth
rate `g`, `cumul(g)` is the long difference of the level.
"""
cumul(y) = CumulResponse(_term(y))

"""
    lags(x, p)
    lags(x, a:b)

Lagged regressors `x_{i,t−1}, …, x_{i,t−p}` (or `x_{i,t−a}, …, x_{i,t−b}`, with `a ≥ 0`).
`x` may be an interaction (`lags(s & X, 1)`) or a first difference
(`lags(firstdiff(y), 1:4)`). Columns are named `L1.x`, `L2.x`, ….
"""
lags(x, p::Integer) = lags(x, 1:Int(p))
function lags(x, r::AbstractUnitRange{<:Integer})
    first(r) >= 0 || throw(ArgumentError("lags must be non-negative; got $r"))
    isempty(r) && throw(ArgumentError("empty lag range $r"))
    return LagsTerm(_term(x), Int(first(r)):Int(last(r)))
end

"""
    firstdiff(x)

First difference `x_{i,t} − x_{i,t−1}` (named `D.x`).
"""
firstdiff(x) = FirstDiffTerm(_term(x))

"""
    fe(v)

Absorbed fixed effect for the groups of column `v`; `fe(a) & fe(b)` absorbs the
interaction (e.g. sector × quarter).
"""
fe(v::Symbol) = FETerm([v])
fe(t::StatsModels.Term) = FETerm([t.sym])

Base.show(io::IO, r::PanelResponse) = print(io, _respname(r))
Base.show(io::IO, t::LagsTerm) = print(io, "lags(", _termlabel(t.term), ", ", t.lags, ")")
Base.show(io::IO, t::FirstDiffTerm) = print(io, "firstdiff(", _termlabel(t.term), ")")
Base.show(io::IO, t::FETerm) = join(io, ("fe($v)" for v in t.vars), " & ")
_termlabel(t) = join(_names(t), ", ")

const _HELPERS = (:leads, :ldiff, :anchor, :cumul, :lags, :firstdiff, :fe)
const _LHS_HELPERS = (:leads, :ldiff, :anchor, :cumul)
# Same-named helpers of these packages are accepted inside `panellp` formulas.
const _FOREIGN_HELPER_MODULES = (:LocalProjections, :Regress, :FixedEffectModels)

# Name of the helper a FunctionTerm stands for, or `nothing`. A same-named function from
# any other module is an error rather than a silently different term.
function _helper(t::StatsModels.FunctionTerm)
    f = t.f
    f isa Function || return nothing
    name = nameof(f)
    name in _HELPERS || return nothing
    owner = nameof(parentmodule(f))
    (owner === :PanelLocalProjections || owner in _FOREIGN_HELPER_MODULES) && return name
    throw(ArgumentError(
        "`$name` in this formula is $(parentmodule(f)).$name, not the panellp helper; " *
        "write PanelLocalProjections.$name(…) in the formula, or rename your function"))
end

# ---------------------------------------------------------------------------------------
# Rewriting
# ---------------------------------------------------------------------------------------

function _lagrange(c::StatsModels.ConstantTerm)
    (p = c.n;
        isinteger(p) && p >= 1 ?
        (1:Int(p)) :
        throw(ArgumentError("lags(x, p) needs an integer p ≥ 1; got $p")))
end
function _lagrange(t::StatsModels.FunctionTerm)
    (t.f === (:) && length(t.args) == 2 &&
     all(a -> a isa StatsModels.ConstantTerm, t.args)) ||
        throw(ArgumentError("lags(x, …) takes an integer p or a range a:b; got $(t.exorig)"))
    a, b = t.args[1].n, t.args[2].n
    (isinteger(a) && isinteger(b) && a >= 0 && b >= a) ||
        throw(ArgumentError("lags(x, a:b) needs integers 0 ≤ a ≤ b; got $(t.exorig)"))
    return Int(a):Int(b)
end
_lagrange(t) = throw(ArgumentError("lags(x, …) takes an integer p or a range a:b"))

# Inner (non-response, non-FE) terms: regressors and the arguments of helpers.
_inner(t::StatsModels.Term) = VarTerm(t.sym)
_inner(t::StatsModels.ContinuousTerm) = VarTerm(t.sym)
_inner(t::StatsModels.ConstantTerm) = ConstTerm(Float64(t.n))
_inner(t::StatsModels.InteractionTerm) = ProductTerm(Any[_inner(c) for c in t.terms])
_inner(t::Union{VarTerm, ConstTerm, FunTerm, ProductTerm, LagsTerm, FirstDiffTerm}) = t
_inner(t::FETerm) = throw(ArgumentError("fe() terms cannot be nested inside other terms"))
function _inner(t::PanelResponse)
    throw(ArgumentError(
        "$(_respname(t)) is a left-hand-side term; horizon-dependent right-hand-side terms are " *
        "not supported"))
end
function _inner(t::StatsModels.FunctionTerm)
    h = _helper(t)
    if h === :lags
        length(t.args) == 2 ||
            throw(ArgumentError("lags takes two arguments, lags(x, p) or lags(x, a:b)"))
        return LagsTerm(_inner(t.args[1]), _lagrange(t.args[2]))
    elseif h === :firstdiff
        length(t.args) == 1 || throw(ArgumentError("firstdiff takes one argument"))
        return FirstDiffTerm(_inner(t.args[1]))
    elseif h === :fe
        throw(ArgumentError("fe() terms cannot be nested inside other terms"))
    elseif h !== nothing
        throw(ArgumentError(
            "$h() is a left-hand-side term; horizon-dependent right-hand-side terms are " *
            "not supported"))
    elseif t.f === (&)
        # `&` inside a helper call is protected by @formula: an interaction of its args
        return ProductTerm(Any[_inner(a) for a in t.args])
    elseif t.f === (~)
        # likewise `~`: an IV block inside a function call
        _inner(StatsModels.FormulaTerm(t.args[1], t.args[2]))
    elseif t.f isa Function && nameof(t.f) in (:lag, :lead) &&
           nameof(parentmodule(t.f)) in (:ShiftedArrays, :StatsModels)
        # ShiftedArrays' lag/lead (re-exported by StatsModels) shift by row position
        throw(ArgumentError(
            "$(nameof(t.f))(x, k) shifts by row position, across units and gaps; in a " *
            "panellp formula write lags(x, k) on the right-hand side (leads are the " *
            "response, leads(y))"))
    end
    return FunTerm(t.f, Any[_inner(a) for a in t.args], string(t.exorig))
end
function _inner(t::StatsModels.FormulaTerm)
    throw(ArgumentError(
        "an instrumental-variable block (endo ~ inst) is a right-hand-side term of its own; " *
        "it cannot be nested inside other terms"))
end
_inner(t) = throw(ArgumentError("unsupported formula term $t"))

function _rewrite_lhs(t)
    t isa PanelResponse && return t
    if t isa StatsModels.FunctionTerm
        h = _helper(t)
        if h === :leads
            length(t.args) == 1 || throw(ArgumentError("leads takes one argument"))
            return LeadsResponse(_inner(t.args[1]))
        elseif h === :ldiff
            length(t.args) == 1 || throw(ArgumentError("ldiff takes one argument"))
            return LDiffResponse(_inner(t.args[1]))
        elseif h === :cumul
            length(t.args) == 1 || throw(ArgumentError("cumul takes one argument"))
            arg = t.args[1]
            if arg isa StatsModels.FunctionTerm && _helper(arg) in _LHS_HELPERS
                throw(ArgumentError("cumul(ldiff(…)) and similar nested responses are not supported"))
            end
            return CumulResponse(_inner(arg))
        elseif h === :anchor
            length(t.args) == 2 || throw(ArgumentError("anchor takes two arguments"))
            return AnchorResponse(_inner(t.args[1]), _inner(t.args[2]))
        elseif h === :fe || h === :lags
            throw(ArgumentError("$h() cannot be the response"))
        elseif t.f === (|)
            throw(ArgumentError(
                "the `|` (anchor pipe) syntax is not supported; use anchor(y, z)"))
        end
    end
    t isa StatsModels.AbstractTerm || throw(ArgumentError("unsupported response $t"))
    t isa StatsModels.InteractionTerm && throw(ArgumentError("unsupported response $t"))
    return PlainResponse(_inner(t))
end

# FE term from `fe(v)` or an interaction of fe() terms; `nothing` if `t` is not one.
function _fe_term(t)
    t isa FETerm && return t
    if t isa StatsModels.FunctionTerm && _helper(t) === :fe
        length(t.args) == 1 && t.args[1] isa StatsModels.Term ||
            throw(ArgumentError("fe() takes one column name, e.g. fe(id)"))
        return FETerm([t.args[1].sym])
    elseif t isa StatsModels.InteractionTerm
        parts = [_fe_term(c) for c in t.terms]
        all(p -> p === nothing, parts) && return nothing
        any(p -> p === nothing, parts) && throw(ArgumentError(
            "fe(a) & x (fixed-effect slopes) is not supported; interact fixed effects only " *
            "with fixed effects, e.g. fe(sector) & fe(quarter)"))
        return FETerm(reduce(vcat, (p.vars for p in parts)))
    end
    return nothing
end

"""
    _rewrite_rhs(rhs) -> (regressors, fes, intercept, iv)

Exogenous regressors in formula order (StatsModels orders main effects before
interactions), the absorbed fixed effects, the intercept request (`:default`, `:yes` for
`1 +`, `:no` for `0 +` / `-1`), and the instrumental-variable block `(endo ~ inst)` as
(endogenous terms, excluded-instrument terms), or `nothing`.
"""
function _rewrite_rhs(rhs)
    terms = rhs isa Tuple ? rhs : (rhs,)
    regs = Any[]
    fes = FETerm[]
    intercept = :default
    iv = nothing
    for t in terms
        if t isa StatsModels.FormulaTerm
            iv === nothing || throw(ArgumentError(
                "the formula has two instrumental-variable blocks; write one, " *
                "(endo1 + endo2 ~ inst1 + inst2 + …), with every endogenous regressor on " *
                "its left and every excluded instrument on its right"))
            iv = _rewrite_iv(t)
            continue
        end
        if t isa StatsModels.ConstantTerm
            if t.n == 1
                intercept = :yes
            elseif t.n == 0 || t.n == -1
                intercept = :no
            else
                throw(ArgumentError("constant $(t.n) is not a valid right-hand-side term"))
            end
            continue
        end
        f = _fe_term(t)
        if f !== nothing
            push!(fes, f)
            continue
        end
        push!(regs, _inner(t))
    end
    return regs, fes, intercept, iv
end

# The block (endo ~ inst) of a panel LP-IV: its endogenous regressors and excluded
# instruments, each side a sum of right-hand-side terms (fixed effects and constants are
# not: they belong to both stages, outside the block).
function _rewrite_iv(t::StatsModels.FormulaTerm)
    function side(x, what)
        out = Any[]
        for c in (x isa Tuple ? x : (x,))
            c isa StatsModels.ConstantTerm && throw(ArgumentError(
                "a constant cannot be $what in (endo ~ inst): the intercept, when there is " *
                "one, is an exogenous regressor of both stages"))
            _fe_term(c) === nothing || throw(ArgumentError(
                "fe() terms cannot be inside (endo ~ inst): absorb fixed effects outside the " *
                "block, where they apply to both stages"))
            push!(out, _inner(c))
        end
        return out
    end
    return (side(t.lhs, "endogenous"), side(t.rhs, "an instrument"))
end

# ---------------------------------------------------------------------------------------
# Names and variables
# ---------------------------------------------------------------------------------------

_names(t::VarTerm) = [string(t.sym)]
_names(t::ConstTerm) = [string(t.value)]
_names(t::FunTerm) = [t.name]
function _names(t::ProductTerm)
    out = _names(t.terms[1])
    for c in t.terms[2:end]
        nc = _names(c)
        out = [string(a, " & ", b) for a in out for b in nc]
    end
    return out
end
_wrap(n::AbstractString) = occursin(' ', n) ? "($n)" : n
function _names(t::LagsTerm)
    [string("L", l, ".", _wrap(n)) for l in t.lags for n in _names(t.term)]
end
_names(t::FirstDiffTerm) = [string("D.", _wrap(n)) for n in _names(t.term)]
_fename(t::FETerm) = join(string.(t.vars), "&")

_respname(r::LeadsResponse) = "leads($(only(_names(r.term))))"
_respname(r::LDiffResponse) = "ldiff($(only(_names(r.term))))"
_respname(r::AnchorResponse) = "anchor($(only(_names(r.term))), $(only(_names(r.base))))"
_respname(r::CumulResponse) = "cumul($(only(_names(r.term))))"
_respname(r::PlainResponse) = only(_names(r.term))

_vars!(s::Set{Symbol}, t::VarTerm) = push!(s, t.sym)
_vars!(s::Set{Symbol}, ::ConstTerm) = s
_vars!(s::Set{Symbol}, t::FunTerm) = (foreach(a -> _vars!(s, a), t.args); s)
_vars!(s::Set{Symbol}, t::ProductTerm) = (foreach(a -> _vars!(s, a), t.terms); s)
_vars!(s::Set{Symbol}, t::Union{LagsTerm, FirstDiffTerm}) = _vars!(s, t.term)
function _vars!(s::Set{Symbol}, r::Union{
        LeadsResponse, LDiffResponse, CumulResponse, PlainResponse})
    _vars!(s, r.term)
end
_vars!(s::Set{Symbol}, r::AnchorResponse) = (_vars!(s, r.term); _vars!(s, r.base))

# Number of response columns must be one.
function _width1(t, what)
    (w = length(_names(t));
        w == 1 ||
            throw(ArgumentError("$what must be a single column; got $w")); t)
end

# ---------------------------------------------------------------------------------------
# Evaluation on the panel: n × w matrices over all data rows, NaN where undefined
# ---------------------------------------------------------------------------------------

_n(ctx::_EvalContext) = nrows(ctx.pidx)

function _cols(t::VarTerm, ctx::_EvalContext)
    haskey(ctx.columns, t.sym) || throw(ArgumentError("column $(t.sym) not found"))
    return reshape(copy(ctx.columns[t.sym]), :, 1)
end
_cols(t::ConstTerm, ctx::_EvalContext) = fill(t.value, _n(ctx), 1)

# A function term is evaluated elementwise on the rows where every data column under it is
# finite, with the natural element types of nested calls (so `ifelse(x > 0, x, 0)` and
# `!(x > 0)` work), and is undefined (NaN) elsewhere: a missing input never becomes a value
# (`NaN > 0` would be a defined `false`).
function _leafcols!(leaves::Vector{Vector{Float64}}, t::FunTerm, ctx::_EvalContext)
    for a in t.args
        a isa ConstTerm && continue
        if a isa FunTerm
            _leafcols!(leaves, a, ctx)
        else
            c = _cols(a, ctx)
            size(c, 2) == 1 ||
                throw(ArgumentError("arguments of $(t.name) must be single columns"))
            push!(leaves, vec(c))
        end
    end
    return leaves
end
function _funeval(t::FunTerm, leaves, rows::Vector{Int}, next::Base.RefValue{Int})
    args = map(t.args) do a
        a isa ConstTerm && return a.value
        a isa FunTerm && return _funeval(a, leaves, rows, next)
        next[] += 1
        return leaves[next[]][rows]
    end
    return t.f.(args...)
end
function _cols(t::FunTerm, ctx::_EvalContext)
    leaves = _leafcols!(Vector{Float64}[], t, ctx)
    rows = [r for r in 1:_n(ctx) if all(c -> isfinite(c[r]), leaves)]
    v = _funeval(t, leaves, rows, Ref(0))
    out = fill(NaN, _n(ctx), 1)
    for (k, r) in enumerate(rows)
        out[r, 1] = Float64(v isa AbstractArray ? v[k] : v)
    end
    return out
end
function _cols(t::ProductTerm, ctx::_EvalContext)
    out = _cols(t.terms[1], ctx)
    for c in t.terms[2:end]
        out = _rowkron(out, _cols(c, ctx))
    end
    return out
end
function _cols(t::LagsTerm, ctx::_EvalContext)
    inner = _cols(t.term, ctx)
    out = Matrix{Float64}(undef, _n(ctx), length(t.lags) * size(inner, 2))
    c = 0
    for l in t.lags
        idx = _shift(ctx.pidx, -l)
        for j in axes(inner, 2)
            c += 1
            out[:, c] = _gather(view(inner, :, j), idx)
        end
    end
    return out
end
function _cols(t::FirstDiffTerm, ctx::_EvalContext)
    inner = _cols(t.term, ctx)
    idx = _shift(ctx.pidx, -1)
    out = similar(inner)
    for j in axes(inner, 2)
        x = view(inner, :, j)
        out[:, j] = x .- _gather(x, idx)
    end
    return out
end

_col1(t, ctx) = vec(_cols(t, ctx))

"""
    _lhs_at(r::PanelResponse, ctx, h) -> Vector{Float64}

The response at horizon `h` for every data row (NaN where undefined).
"""
function _lhs_at(r::LeadsResponse, ctx::_EvalContext, h::Integer)
    _gather(_col1(r.term, ctx), _shift(ctx.pidx, h))
end
function _lhs_at(r::LDiffResponse, ctx::_EvalContext, h::Integer)
    y = _col1(r.term, ctx)
    return _gather(y, _shift(ctx.pidx, h)) .- _gather(y, _shift(ctx.pidx, -1))
end
function _lhs_at(r::AnchorResponse, ctx::_EvalContext, h::Integer)
    _gather(_col1(r.term, ctx), _shift(ctx.pidx, h)) .- _col1(r.base, ctx)
end
function _lhs_at(r::CumulResponse, ctx::_EvalContext, h::Integer)
    y = _col1(r.term, ctx)
    acc = _gather(y, _shift(ctx.pidx, 0))
    for k in 1:h
        acc .+= _gather(y, _shift(ctx.pidx, k))
    end
    return acc
end
function _lhs_at(r::PlainResponse, ctx::_EvalContext, h::Integer)
    h == 0 || throw(ArgumentError(
        "a plain response is a pre-built outcome and is estimated at horizon 0 only; use " *
        "leads(y), ldiff(y), anchor(y, z) or cumul(y) to estimate several horizons"))
    return _col1(r.term, ctx)
end
