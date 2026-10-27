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

TLAHR(; refine::Bool = true, weights::Symbol = :residual) = TLAHR(refine, weights)

function CosineHAR(B = :llsw)
    B isa Integer && return CosineHAR(:fixed, Int(B))
    B === :llsw || throw(ArgumentError("B must be :llsw or a number of cosine terms"))
    return CosineHAR(:llsw, 0)
end

# Number of cosine terms of CosineHAR(:llsw) for T periods: ⌊0.41 T^{2/3}⌋, as harreg
# (har_utils.mata), computed the same way in double precision (T = 1000 gives 40).
_llsw(T::Int) = floor(Int, 0.41 * T^(2 / 3))

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
function _label(s::TLAHR)
    opts = String[]
    s.refine || push!(opts, "refine = false")
    s.weights === :residual || push!(opts, "weights = :$(s.weights)")
    return "TLAHR($(join(opts, ", ")))"
end
_label(s::CosineHAR) = "CosineHAR($(s.rule === :llsw ? ":llsw" : s.B))"

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

# t-LAHR without the refinement: CR0 by period, normal critical values (A&S eq. 11);
# defined where Clustered(:time; ssc = :none) is.
function _vcov_parts(s::TLAHR, r, fit::HorizonFit, ctx)
    T = count(fit.timepresent)
    defined = T > 1 && fit.dof.nobs - fit.dof.k - fit.dof.df_a > 0
    c = _finite_or_nan(1.0, defined)
    return _meat1(fit, 2), c, defined ? Inf : NaN, VcovInfo(c, T, T, 0, false)
end

# Number of cosine terms and calendar periods (first to last of the sample) of a fit.
function _cosineterms(s::CosineHAR, fit::HorizonFit)
    T = findlast(fit.timepresent) - findfirst(fit.timepresent) + 1
    return s.rule === :llsw ? _llsw(T) : s.B, T
end

# Fewer cosine terms than periods are needed (B = ⌊0.41 T^{2/3}⌋ is 0 for T ≤ 3).
function _sandwich(s::CosineHAR, r, fit::HorizonFit, ctx, kidx::Vector{Int})
    B, T = _cosineterms(s, fit)
    nk = length(kidx)
    1 <= B < T ||
        return fill(NaN, nk, nk), fill(NaN, nk), VcovInfo(NaN, 0, T, B, false), :bandwidth
    return _plain_sandwich(s, r, fit, ctx, kidx)
end

# Equal-weighted cosine variance of the period sums of the scores, over the calendar
# periods from the first to the last of the sample (harreg's EWC: Λ_j = √2 Σ_t
# cos(πj(t − ½)/T) S_t, Ω = Σ_j Λ_j Λ_j' / B, no small-sample factor, t(B)).
function _vcov_parts(s::CosineHAR, r, fit::HorizonFit, ::Nothing)
    present = fit.timepresent
    g = fit.scores.S[2][findfirst(present):findlast(present), :]
    B, T = _cosineterms(s, fit)
    defined = 1 <= B < T && fit.dof.nobs - fit.dof.k - fit.dof.df_a > 0
    M = zeros(size(g, 2), size(g, 2))
    if defined
        Φ = [sqrt(2) * cos(π * j * (t - 0.5) / T) for j in 1:B, t in 1:T]
        Λ = Φ * g
        M = (Λ' * Λ) ./ B
    end
    c = _finite_or_nan(1.0, defined)
    return M, c, defined ? Float64(B) : NaN, VcovInfo(c, 0, T, B, false)
end

# t-LAHR: every coefficient uses the synthetic regression of its owner (`r.owner`, from
# `_owners`). With the refinement, every estimated coefficient needs one, so every shock
# must be an aggregate series times an exposure.
function _prepare(s::TLAHR, r::PanelLocalProjection)
    r.estimator isa SplitPanelJackknife && throw(ArgumentError(
        "TLAHR is the inference of a lag-augmented fixed-effects regression, not of the " *
        "split-panel jackknife; use Clustered or DriscollKraaySE"))
    s.refine && !isempty(r.endogenous) &&
        throw(ArgumentError(
            "TLAHR(refine = true) is not available for panel LP-IV yet (Almuzara & Sancibrián " *
            "define the refinement for OLS only); TLAHR(refine = false) is their LP-IV " *
            "variance, CR0 by period on the IV scores (supplement A.4)"))
    if s.refine
        for f in r.fits
            f.status === :ok || continue
            for j in findall(f.keep)
                o = r.owner[j]
                (o > 0 && any(sy -> sy.col == o, f.synthetic)) && continue
                bad = [nm
                       for (k, nm) in enumerate(r.coefnames)
                       if nm in r.shocks && r.owner[k] == 0]
                which = join(repr.(isempty(bad) ? r.shocks[1:1] : bad), ", ")
                throw(ArgumentError(
                    "TLAHR(refine = true) needs every shock to be an aggregate " *
                    "series, or an aggregate series times a unit-level exposure " *
                    "(s & X); $which is neither. TLAHR(refine = false) is the CR0 " *
                    "variance clustered by period"))
            end
        end
    end
    return (; owner = r.owner)
end

"""
    _ikdesign(X) -> NamedTuple or nothing

Synthetic regression of the HC2/Imbens–Kolesár refinement (A&S supplement B.2): from the
thin SVD of the synthetic regressors (T × K; their cross-product is never formed), keeping
the singular values above max(T, K)·eps·σ₁, the leverages h_t and the residual-maker
diagonal d_t = 1 − h_t. `nothing` when a period has leverage one (d_t ≈ 0) or a regressor
is not finite.
"""
function _ikdesign(X::Matrix{Float64})
    all(isfinite, X) || return nothing
    T, K = size(X)
    F = LinearAlgebra.svd(X)
    tol = max(T, K) * eps() * (isempty(F.S) ? 0.0 : F.S[1])
    r = count(>(tol), F.S)
    U = F.U[:, 1:r]
    h = vec(sum(abs2, U; dims = 2))
    d = 1 .- h
    minimum(d; init = 1.0) > 1e-10 || return nothing
    return (; d, h, U, Sr = F.S[1:r], V = F.V[:, 1:r])
end

# Bell–McCaffrey/Imbens–Kolesár degrees of freedom of coefficient c (column of X̄), with the
# homoskedastic reference: ν = (Σλ)²/Σλ² for the eigenvalues λ of G'G, G_{·t} = M̄_{·t}
# a_t/√d_t, a = X̄(X̄'X̄)⁺e_c, M̄ = I − UU'. Since G'G = D M̄ D with D = diag(a/√d),
# Σλ = Σ_t a_t² and Σλ² = Σ_{t,s} b_t b_s M̄_ts², b = a²/d: the diagonal Σ_t b_t² d_t², and
# off it Σ_{t≠s} b_t b_s (UU')_ts² = ‖U' diag(b) U‖²_F − Σ_t b_t² h_t². No
# eigendecomposition, and no T × T matrix. The contrast is e_c of each coefficient (the
# authors' code uses e₁ for every coefficient).
function _iknu(des, c::Int)
    a = des.U * (des.V[c, :] ./ des.Sr)
    b = a .^ 2 ./ des.d
    off = sum(abs2, des.U' * (b .* des.U)) - sum(abs2, b .* des.h)
    den = sum(abs2, b .* des.d) + max(off, 0.0)
    return den > 0 ? sum(abs2, a)^2 / den : NaN
end

# Variance, per-coefficient critical df, details and flag (:ok, :undefined, :negfactor or
# :leverage) of one horizon, on the kept coefficients `kidx`.
function _sandwich(spec::PanelVcov, r, fit::HorizonFit, ctx, kidx::Vector{Int})
    _plain_sandwich(spec, r, fit, ctx, kidx)
end

# c · B M B with the meat, factor and df of `_vcov_parts`.
function _plain_sandwich(spec::PanelVcov, r, fit::HorizonFit, ctx, kidx::Vector{Int})
    M, c, df, info = _vcov_parts(spec, r, fit, ctx)
    nk = length(kidx)
    isnan(c) && return fill(NaN, nk, nk), fill(NaN, nk), info, :undefined
    B = fit.bread[kidx, kidx]
    Vk = c .* (B * M[kidx, kidx] * B)
    return (Vk + Vk') ./ 2, fill(df, nk), info, c < 0 ? :negfactor : :ok
end

"""
t-LAHR with the refinement (A&S supplement B.2). Each owner k (`_owners`: a shock, or
another aggregate series times an exposure) has its synthetic regression (stored by the
fit): the per-period means X̄⁽ᵏ⁾_t weighted by ω⁽ᵏ⁾_t (fn 24), so with rows √ω_t X̄_t, whose
leverage is that of the period in the panel regression when every regressor is an
aggregate series times the exposure (HC2 is then CR2 clustered by period); unweighted X̄_t
with the raw exposure for `weights = :raw`, as the authors' code. With residual-maker
diagonal d⁽ᵏ⁾_t and z⁽ʲ⁾_t = (B S_t)_j / √d⁽ᶜ⁽ʲ⁾⁾_t for coefficient j of owner c(j),
V_jl = Σ_t z⁽ʲ⁾_t z⁽ˡ⁾_t: the paper's HC2 variance B (Σ_t S_t S_t'/d_t) B for every
coefficient's own variance, and positive semi-definite. Each coefficient has its own
Imbens–Kolesár df, from `_iknu`.
"""
function _sandwich(s::TLAHR, r, fit::HorizonFit, ctx::NamedTuple, kidx::Vector{Int})
    s.refine || return _plain_sandwich(s, r, fit, ctx, kidx)
    T = count(fit.timepresent)
    nk = length(kidx)
    info = VcovInfo(1.0, T, T, 0, false)
    nan = (fill(NaN, nk, nk), fill(NaN, nk), VcovInfo(NaN, T, T, 0, false))
    (T > 1 && fit.dof.nobs - fit.dof.k - fit.dof.df_a > 0) || return (nan..., :undefined)
    S = fit.scores.S[2][fit.timepresent, kidx]
    G = S * fit.bread[kidx, kidx]
    designs = Dict{Int, Any}()
    for o in unique(ctx.owner[kidx])
        sy = fit.synthetic[findfirst(x -> x.col == o, fit.synthetic)]
        X̄ = s.weights === :raw ? sy.Xraw[:, kidx] : sqrt.(sy.weight) .* sy.X[:, kidx]
        des = _ikdesign(X̄)
        des === nothing && return (nan..., :leverage)
        designs[o] = des
    end
    Z = similar(G)
    dfs = fill(NaN, nk)
    for (c, j) in enumerate(kidx)
        des = designs[ctx.owner[j]]
        Z[:, c] = G[:, c] ./ sqrt.(des.d)
        dfs[c] = _iknu(des, c)
    end
    return Z' * Z, dfs, info, :ok
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
[`Clustered`](@ref), [`DriscollKraaySE`](@ref), [`TLAHR`](@ref) or [`CosineHAR`](@ref)
specification. Omitted coefficients have NaN rows and columns. Horizons where the variance
is undefined (a single cluster or period, no residual degrees of freedom, or a period with
leverage one in the t-LAHR synthetic regression) are NaN, with a warning; so are, in the
standard errors, negative variances of a two-way variance without PSD repair.
"""
function StatsAPI.vcov(spec::PanelVcov, r::PanelLocalProjection)
    K = length(r.coefnames)
    ctx = _prepare(spec, r)
    metric = _psdmetric(spec)
    Vs = Matrix{Float64}[]
    dfs = Vector{Float64}[]
    infos = VcovInfo[]
    undefined, leverage, negfactor, negvar = Int[], Int[], Int[], Int[]
    bandwidth = Tuple{Int, Int, Int}[]
    for fit in r.fits
        V = fill(NaN, K, K)
        if fit.status !== :ok
            push!(Vs, V)
            push!(dfs, fill(NaN, K))
            push!(infos, VcovInfo(NaN, 0, 0, 0, false))
            continue
        end
        kidx = findall(fit.keep)
        Vk, dk, info, flag = _sandwich(spec, r, fit, ctx, kidx)
        fixed = false
        if flag === :bandwidth
            push!(bandwidth, (fit.h, info.lags, info.periods))
        elseif flag === :undefined || flag === :leverage
            push!(flag === :undefined ? undefined : leverage, fit.h)
        else
            flag === :negfactor && push!(negfactor, fit.h)
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
        d[kidx] .= dk
        push!(dfs, d)
        push!(infos, VcovInfo(info.factor, info.groups, info.periods, info.lags, fixed))
    end
    label = _label(spec)
    isempty(bandwidth) ||
        @warn "$label needs fewer cosine terms than periods (1 ≤ B < T, T the periods " *
              "from the first to the last of the sample): B = $(bandwidth[1][2]) and " *
              "T = $(bandwidth[1][3]) at h = $(bandwidth[1][1]); it is NaN at " *
              "h = $(_hlist(first.(bandwidth)))"
    isempty(undefined) ||
        @warn "$label is undefined at h = $(_hlist(undefined)) (a single " *
              "cluster or period, or no residual degrees of freedom); it is NaN there"
    isempty(leverage) ||
        @warn "$label: a period has leverage one in the synthetic regression at " *
              "h = $(_hlist(leverage)) (too few periods for the regressors?); the " *
              "variance is NaN there"
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

"""
    firststage(r; vcov = TLAHR(refine = false)) -> FirstStageSummary

First-stage diagnostics of a panel LP-IV result, per horizon, with the variance
specification `vcov` (any of [`Robust`](@ref), [`Clustered`](@ref),
[`DriscollKraaySE`](@ref), [`CosineHAR`](@ref), `TLAHR(refine = false)`) applied to the
first-stage scores, with its small-sample factor computed from the first stage's degrees
of freedom (its regressors: the exogenous ones and the excluded instruments), as ivreg2,
ivreghdfe and weakivtest do:

- `kpF`: the Kleibergen–Paap (2006) rk Wald F for weak identification, rk/L, with rk the
  Wald statistic of the test that Π (instruments × endogenous, net of the fixed effects and
  the exogenous regressors) has rank m − 1, from the singular value decomposition of
  Θ̂ = G Π̂ F′, G'G = Z̃'Z̃/n, F'F = (X̃'X̃/n)⁻¹ (their eqs. 11–14). With one endogenous
  regressor it is the Wald F of Π = 0. ivreg2's `widstat`.
- `effF` (one endogenous regressor): the Montiel Olea & Pflueger (2013) effective F,
  π̂'Z̃'Z̃π̂ / tr(V̂(π̂) Z̃'Z̃); equal to `kpF` with one instrument.
- `coef`, `se`: the first-stage coefficients of the instruments and their standard errors
  (instruments × endogenous) per horizon.

Two-way clustered variances are not repaired (`psd` is ignored, as in ivreg2): where the
variance that a statistic needs is not positive definite, the statistic is NaN, with a
warning.
"""
function firststage(r::PanelLocalProjection; vcov::PanelVcov = TLAHR(refine = false))
    isempty(r.endogenous) &&
        throw(ArgumentError("firststage needs a panel LP-IV result (an (endo ~ inst) block)"))
    ctx = _prepare(vcov, r)
    H = length(r.fits)
    kpF, effF, rk = fill(NaN, H), fill(NaN, H), fill(NaN, H)
    ninst, nendo, nobs = zeros(Int, H), zeros(Int, H), zeros(Int, H)
    endos, insts = [String[] for _ in 1:H], [String[] for _ in 1:H]
    coefs, ses = [zeros(0, 0) for _ in 1:H], [zeros(0, 0) for _ in 1:H]
    undefined = Int[]
    for (k, fit) in enumerate(r.fits)
        nobs[k] = fit.dof.nobs
        fs = fit.firststage
        (fit.status === :ok && fs !== nothing) || continue
        L, m = size(fs.coef)
        ninst[k], nendo[k] = L, m
        endos[k] = r.coefnames[fs.endo]
        insts[k] = r.instruments[fs.inst]
        coefs[k] = fs.coef
        Lm = L * m
        # the first stage's regressors: the structural ones less the endogenous, plus the
        # excluded instruments
        dv = fit.dof
        dof = DofInfo(dv.nobs, dv.k - m + L, dv.df_a, dv.felevels, dv.fecomponents,
            dv.ngroups, dv.nested, dv.singletons)
        pfit = HorizonFit(fit.h, :ok, vec(fs.coef), trues(Lm), falses(Lm), fs.bread,
            ones(Lm), fit.esample, fit.timepresent, dof, fs.scores, fit.iterations,
            fit.converged, trues(Lm), fit.nlags, SyntheticSeries[], falses(0), nothing)
        V, _, _, flag = _sandwich(vcov, r, pfit, ctx, collect(1:Lm))
        if flag === :undefined || flag === :bandwidth || !all(isfinite, V)
            push!(undefined, fit.h)
            continue
        end
        ses[k] = reshape([v >= 0 ? sqrt(v) : NaN for v in LinearAlgebra.diag(V)], L, m)
        rk[k] = _kprk(fs, V, fit.dof.nobs)
        isnan(rk[k]) && push!(undefined, fit.h)
        kpF[k] = rk[k] / L
        if m == 1
            π = vec(fs.coef)
            den = LinearAlgebra.tr(V * fs.zz)
            effF[k] = den > 0 ? (π' * fs.zz * π) / den : NaN
        end
    end
    isempty(undefined) || @warn "firststage: $(_label(vcov)) is undefined, or the " *
          "variance of the first stage is not positive definite, at " *
          "h = $(_hlist(sort!(unique(undefined)))); the statistics are NaN there"
    return FirstStageSummary(_label(vcov), collect(r.horizons), kpF, effF, rk, ninst,
        nendo, nobs, endos, insts, coefs, ses)
end

# Kleibergen–Paap (2006) rk Wald statistic of the test that Π (L × m) has rank m − 1, with
# V = Var(vec Π̂): Θ̂ = L_z' Π̂ L_y⁻ᵀ (L_z L_z' = Z̃'Z̃/n, L_y L_y' = X̃'X̃/n), Θ̂ = U S V' (full
# SVD), q = m − 1; λ = (B ⊗ A') vec Θ̂ and rk = λ' [(B ⊗ A') V_Θ (B ⊗ A')']⁻¹ λ, with
# V_Θ = (L_y⁻¹ ⊗ L_z') V (L_y⁻¹ ⊗ L_z')'. KP's eq. (12) takes A = U₂ U₂₂⁻¹(U₂₂U₂₂')^{1/2}
# and B = (V₂₂V₂₂')^{1/2} V₂₂'⁻¹ V₂' (U₂, V₂ the trailing singular vectors): nonsingular
# transformations of U₂ and V₂, which leave the Wald statistic unchanged. With A = U₂ and
# B = V₂' there is nothing to invert, also when U₂₂ or V₂₂ is singular (tied singular
# values, e.g. group dummies × an aggregate with a diagonal Π̂).
function _kprk(fs::FirstStage, V::Matrix{Float64}, n::Int)
    L, m = size(fs.coef)
    Lz = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(fs.zz ./ n)).L
    Ly = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(fs.yy ./ n)).L
    iLy = inv(Ly)
    Θ = Lz' * fs.coef * iLy'
    F = LinearAlgebra.svd(Θ; full = true)
    A = F.U[:, m:L]                  # trailing left singular vectors (q + 1 = m onward)
    B = F.V[:, m:m]'                 # the last right singular vector
    T = kron(iLy, Lz')
    K = kron(B, A')
    λ = K * vec(Θ)
    Vλ = LinearAlgebra.Symmetric(K * (T * V * T') * K')
    # the Wald statistic needs Var(λ̂) positive definite; an indefinite one (two-way
    # clustering is not repaired) gives NaN, never a negative statistic
    ev = LinearAlgebra.eigvals(Vλ)
    first(ev) > 1e-12 * max(last(ev), 0.0) || return NaN
    return only(λ' * (Vλ \ λ))
end
