# Variance estimators, computed after estimation from the score sufficient statistics of
# each horizon (scores.jl). Inspired by CovarianceMatrices.jl (G. Ragusa) and the
# variance code of reghdfe (S. Correia) and xtlp (Mei, Sheng & Shi); written for the
# stored group sums. Two departures from CovarianceMatrices on purpose: the
# Driscoll–Kraay sums are ordered by calendar period (not first appearance), and the CR1
# factor multiplies the variance once (not the scores, which would square G/(G−1)).

# Keyword constructors; the inner constructors in types.jl validate every field.
Robust(; ssc::Symbol = :stata) = Robust(ssc)

function Clustered(dims...; ssc::Symbol = :stata, psd = length(dims) == 2)
    1 <= length(dims) <= 2 ||
        throw(ArgumentError("Clustered takes one or two dimensions; got $(length(dims))"))
    names = map(dims) do d
        d isa Symbol && return d
        d isa Tuple{Vararg{Symbol}} && return _clustername(d)
        throw(ArgumentError("cluster dimensions are :unit, :time, column names or tuples of " *
                            "column names; got $d"))
    end
    p = psd === true ? :stata : psd === false ? :none : psd
    p isa Symbol || throw(ArgumentError("psd must be true, false, :stata, :raw or :none"))
    return Clustered(collect(Symbol, names), ssc, p)
end

function DriscollKraaySE(bw = :nw94; ssc::Symbol = :stata, timelags::Symbol = :calendar)
    bw isa Integer && return DriscollKraaySE(:fixed, Int(bw), ssc, timelags)
    (bw isa Symbol && bw !== :fixed) ||
        throw(ArgumentError("bw must be :nw94, :as2026, :mss2026 or a number of lags"))
    return DriscollKraaySE(bw, 0, ssc, timelags)
end

"""
    _dk_lags(rule, T) -> m

Driscoll–Kraay lags for `T` time periods. `:as2026` is Almuzara–Sancibrián's kernel
bandwidth M = ⌈4(T/100)^{2/9}⌉ with Bartlett weights 1 − ℓ/M, i.e. M − 1 lags with
weights 1 − j/(m+1).
"""
function _dk_lags(s::DriscollKraaySE, T::Int)
    x = 4 * (T / 100)^(2 / 9)
    m = s.rule === :nw94 ? floor(Int, x) :
        s.rule === :as2026 ? ceil(Int, x) - 1 :
        s.rule === :mss2026 ? floor(Int, T^(1 / 4)) : s.lags
    return clamp(m, 0, max(T - 1, 0))
end

_label(s::Robust) = "Robust(ssc = :$(s.ssc))"
function _label(s::Clustered)
    "Clustered($(join(repr.(s.dims), ", ")); ssc = :$(s.ssc)$(s.psd === :none || length(s.dims) == 1 ? "" : ", psd = :$(s.psd)"))"
end
function _label(s::DriscollKraaySE)
    "DriscollKraaySE($(s.rule === :fixed ? s.lags : repr(s.rule)); ssc = :$(s.ssc)$(s.timelags === :position ? ", timelags = :position" : ""))"
end

function _dimindex(r::PanelLocalProjection, d::Symbol)
    nm = d === r.panel.id ? :unit : d === r.panel.time ? :time : d
    g = findfirst(s -> _samedim(nm, s), r.groupnames)
    if g === nothing
        parts = split(string(d), '&')
        hint = length(parts) == 1 ? "clusters = (:$d,)" :
               "clusters = (($(join(":" .* parts, ", "))),)"
        throw(ArgumentError(
            "no score sums stored for cluster dimension :$d. Stored: " *
            "$(join(repr.(r.groupnames), ", ")). Pass $hint to panellp"))
    end
    return g
end

# Score dimensions of a specification, resolved once per `vcov` call.
_prepare(::PanelVcov, r) = nothing
function _prepare(s::Clustered, r::PanelLocalProjection)
    gs = [_dimindex(r, d) for d in s.dims]
    length(gs) == 2 && gs[1] == gs[2] &&
        throw(ArgumentError(
            "the two cluster dimensions $(repr(s.dims[1])) and $(repr(s.dims[2])) are the same"))
    if s.ssc === :mss
        names = Set(r.groupnames[g] for g in gs)
        names == Set([:unit]) || names == Set([:unit, :time]) ||
            throw(ArgumentError(
                "ssc = :mss is the Mei–Sheng–Shi factor for Clustered(:unit) and " *
                "Clustered(:unit, :time) only; use ssc = :stata or :none for $(_label(s))"))
    end
    return gs
end

_meat1(fit::HorizonFit, g::Int) = (S = fit.scores.S[g]; S' * S)

function _meat2(fit::HorizonFit, a::Int, b::Int)
    a > b && ((a, b) = (b, a))
    st = fit.scores
    k = findfirst(==((a, b)), st.pairs)
    kind = st.pairkind[k]
    Ma = _meat1(fit, a)
    Mb = _meat1(fit, b)
    kind === :first && return Mb      # a ⊆ b: M_a + M_b − M_a
    kind === :second && return Ma
    Mab = kind === :obs ? st.Mobs : (Sab = st.Spair[k]; Sab' * Sab)
    return Ma + Mb - Mab
end

# Absorbed degrees of freedom when the fixed effects nested in any of the score dimensions
# `gs` are dropped from the count (reghdfe's treatment of FE nested in a cluster
# variable); returns (df_a of the remaining FE, whether any was nested).
function _dfa_nested(fit::HorizonFit, gs::Vector{Int})
    D = size(fit.dof.nested, 1)
    nested = [any(g -> fit.dof.nested[f, g], gs) for f in 1:D]
    remaining = [f for f in 1:D if !nested[f]]
    anyn = any(nested)
    return _df_fe(fit.dof.felevels, fit.dof.fecomponents, remaining, anyn), anyn
end

# (meat, factor, critical df, info) of one horizon. The factor is NaN where the variance
# is undefined: a single cluster or period, or no residual degrees of freedom.
_finite_or_nan(c, ok::Bool) = ok ? Float64(c) : NaN

function _vcov_parts(s::Robust, r, fit::HorizonFit, ::Nothing)
    n, k, dfa = fit.dof.nobs, fit.dof.k, fit.dof.df_a
    dfr = n - k - dfa
    defined = dfr > 0
    c, df = s.ssc === :none ? (1.0, Inf) : (n / dfr, Float64(dfr))
    c = _finite_or_nan(c, defined)
    return fit.scores.Mobs, c, defined ? df : NaN, VcovInfo(c, 0, 0, 0, false)
end

function _vcov_parts(s::Clustered, r, fit::HorizonFit, gs::Vector{Int})
    M = length(gs) == 1 ? _meat1(fit, gs[1]) : _meat2(fit, gs[1], gs[2])
    n, k, dfa = fit.dof.nobs, fit.dof.k, fit.dof.df_a
    G = minimum(fit.dof.ngroups[g] for g in gs)
    defined = G > 1 && n - k - dfa > 0
    if s.ssc === :none
        c, df = 1.0, Inf
    elseif s.ssc === :stata
        # reghdfe (and xtlp for one or two FE): (n−1)/(n − a − K − df_a_rem) · G/(G−1),
        # df = min(n − K − df_a_rem, G − 1)
        dfr, anyn = _dfa_nested(fit, gs)
        defined &= n - anyn - k - dfr > 0
        c = (n - 1) / (n - anyn - k - dfr) * G / (G - 1)
        df = Float64(min(n - k - dfr, G - 1))
    else  # :mss, Mei–Sheng–Shi replication code (LP_panel_all.r, robust = T)
        N, T0 = r.panel.nunits, r.panel.ntimes
        if length(gs) == 1
            c = N / (N - 1) * (n - 1) / (n - k)
            G = N
        else
            # N counts every unit in the data: the factor turns negative when n ≤ N + K
            # (reported by vcov), as in the authors' code
            G = min(N, T0)
            c = G / (G - 1) * (n - 1) / (n - N - k)
        end
        defined &= isfinite(c)
        df = Inf
    end
    c = _finite_or_nan(c, defined)
    return M, c, defined ? df : NaN, VcovInfo(c, G, 0, 0, false)
end

function _vcov_parts(s::DriscollKraaySE, r, fit::HorizonFit, ::Nothing)
    St = fit.scores.S[2]                      # time sums, calendar order
    present = fit.timepresent
    g = s.timelags === :calendar ? St : St[present, :]
    T = count(present)
    m = _dk_lags(s, T)
    M = g' * g
    for j in 1:m
        Γ = view(g, (j + 1):size(g, 1), :)' * view(g, 1:(size(g, 1) - j), :)
        M .+= (1 - j / (m + 1)) .* (Γ .+ Γ')
    end
    n, k, dfa = fit.dof.nobs, fit.dof.k, fit.dof.df_a
    defined = T > 1 && n - k - dfa > 0
    if s.ssc === :none
        c, df = 1.0, Inf
    elseif s.ssc === :stata
        # reghdfe/xtlp: time is the cluster variable for the dof count, without nested_adj
        dfdk, _ = _dfa_nested(fit, [2])
        defined &= n - k - dfdk > 0
        c = (n - 1) / (n - k - dfdk) * T / (T - 1)
        df = Float64(min(n - k - dfdk, T - 1))
    else  # :xtscc
        c = (n - 1) / (n - k - (dfa > 0 ? 1 : 0)) * T / (T - 1)
        df = Float64(T - 1)
    end
    c = _finite_or_nan(c, defined)
    return M, c, defined ? df : NaN, VcovInfo(c, 0, T, m, false)
end

_psdmetric(s::Clustered) = length(s.dims) == 2 ? s.psd : :none
_psdmetric(s::PanelVcov) = :none

# Clip negative eigenvalues of V to zero (Cameron, Gelbach & Miller 2011). With `sd`, the
# clipping is done on S V S, S = diag(sd), and mapped back (reghdfe/xtlp metric). With an
# intercept at position `icol`, the slope block is then replaced by its own clip, as
# reghdfe does when it reports `_cons` (reghdfe.mata, `reghdfe_fix_psd`).
function _psd_clip(V::Matrix{Float64}, sd::Union{Nothing, Vector{Float64}},
        icol::Union{Nothing, Int} = nothing)
    W = sd === nothing ? V : V .* (sd * sd')
    Wc, fixed = _clip(W)
    fixed || return V, false
    if icol !== nothing
        s = [j for j in axes(W, 1) if j != icol]
        Wc[s, s] = _clip(W[s, s])[1]
    end
    return (sd === nothing ? Wc : Wc ./ (sd * sd')), true
end
function _clip(W::Matrix{Float64})
    E = LinearAlgebra.eigen(LinearAlgebra.Symmetric(W))
    any(<(0), E.values) || return W, false
    return E.vectors * LinearAlgebra.Diagonal(max.(E.values, 0.0)) * E.vectors', true
end

"""
    vcov(spec::PanelVcov, r::PanelLocalProjection) -> PanelLPVcov

Variance of the coefficients at every horizon, for a [`Robust`](@ref),
[`Clustered`](@ref) or [`DriscollKraaySE`](@ref) specification. Omitted coefficients have
NaN rows and columns. Horizons where the variance is undefined (a single cluster or
period, or no residual degrees of freedom) are NaN, with a warning; so are, in the
standard errors, negative variances of a two-way variance without PSD repair.
"""
function StatsAPI.vcov(spec::PanelVcov, r::PanelLocalProjection)
    K = length(r.coefnames)
    ctx = _prepare(spec, r)
    metric = _psdmetric(spec)
    Vs = Matrix{Float64}[]
    dfs = Vector{Float64}[]
    infos = VcovInfo[]
    undefined, negfactor, negvar = Int[], Int[], Int[]
    for fit in r.fits
        V = fill(NaN, K, K)
        if fit.status !== :ok
            push!(Vs, V)
            push!(dfs, fill(NaN, K))
            push!(infos, VcovInfo(NaN, 0, 0, 0, false))
            continue
        end
        kidx = findall(fit.keep)
        M, c, df, info = _vcov_parts(spec, r, fit, ctx)
        fixed = false
        if isnan(c)
            push!(undefined, fit.h)
        else
            c < 0 && push!(negfactor, fit.h)
            B = fit.bread[kidx, kidx]
            Vk = c .* (B * M[kidx, kidx] * B)
            Vk = (Vk + Vk') ./ 2
            if metric !== :none
                icol = r.intercept && fit.keep[1] ? 1 : nothing
                Vk, fixed = metric === :stata ? _psd_clip(Vk, fit.xsd[kidx], icol) :
                            _psd_clip(Vk, nothing)
            end
            any(<(0), LinearAlgebra.diag(Vk)) && push!(negvar, fit.h)
            V[kidx, kidx] = Vk
        end
        push!(Vs, V)
        d = fill(NaN, K)
        d[kidx] .= df
        push!(dfs, d)
        push!(infos, VcovInfo(info.factor, info.groups, info.periods, info.lags, fixed))
    end
    label = _label(spec)
    isempty(undefined) ||
        @warn "$label is undefined at h = $(_hlist(undefined)) (a single " *
              "cluster or period, or no residual degrees of freedom); it is NaN there"
    isempty(negfactor) || @warn "$label: the Mei–Sheng–Shi two-way factor is negative at " *
          "h = $(_hlist(negfactor)) (n ≤ N + K, with N the units in the data)"
    isempty(negvar) || @warn "$label has negative variances at h = $(_hlist(negvar)); " *
          "their standard errors are NaN. A PSD repair (psd = :stata or :raw) avoids them"
    return PanelLPVcov(spec, r.coefnames, r.horizons, Vs, dfs, infos, _fingerprint(r))
end

StatsAPI.vcov(cv::PanelLPVcov, h::Integer) = cv.V[_hindex(cv.horizons, h)]

function _hindex(horizons::UnitRange{Int}, h::Integer)
    h in horizons ||
        throw(ArgumentError("horizon $h not estimated; horizons are $horizons"))
    return h - first(horizons) + 1
end

# Identifies the result a PanelLPVcov was computed from.
function _fingerprint(r::PanelLocalProjection)
    h = hash(r.coefnames, hash(r.horizons))
    for f in r.fits
        h = hash(f.dof.nobs, hash(f.status, h))
        for x in f.coef
            h = hash(x, h)
        end
    end
    return h
end
