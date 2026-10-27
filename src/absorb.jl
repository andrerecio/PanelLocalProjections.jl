# Fixed-effect absorption on one horizon's estimation sample.
#
# Design follows the within/alternating-projections approach of FixedEffectModels.jl
# (M. Gomez) and Regress.jl (G. Ragusa; src/utils/fit_common.jl, utils/fixedeffects.jl)
# and the Irons–Tuck acceleration of fixest (L. Bergé) and reghdfe (S. Correia); no code is
# copied from them. The degrees-of-freedom rule for two dimensions (levels minus connected
# components of the bipartite group graph) is the one reghdfe uses for mobility groups.

"""
    _drop_singletons!(mask, fe) -> Int

Remove from `mask` the rows that are alone in their group of some fixed-effect dimension,
repeatedly, until no singleton is left (the order of removal does not change the result).
Returns the number of rows removed.
"""
function _drop_singletons!(mask::BitVector, fe::Vector{Vector{Int32}})
    isempty(fe) && return 0
    counts = [zeros(Int, maximum(c; init = Int32(0))) for c in fe]
    dropped = 0
    while true
        for (d, c) in enumerate(fe)
            cnt = counts[d]
            fill!(cnt, 0)
            @inbounds for r in eachindex(mask)
                mask[r] && (cnt[c[r]] += 1)
            end
        end
        changed = false
        @inbounds for r in eachindex(mask)
            mask[r] || continue
            for (d, c) in enumerate(fe)
                if counts[d][c[r]] == 1
                    mask[r] = false
                    dropped += 1
                    changed = true
                    break
                end
            end
        end
        changed || return dropped
    end
end

"""
    _fe_sample(fe, rows) -> _FESample

Fixed-effect codes on the estimation rows, re-indexed to `1:L` per dimension in order of
first appearance, with group sizes.
"""
function _fe_sample(fe::Vector{Vector{Int32}}, rows::Vector{Int})
    codes = Vector{Vector{Int32}}(undef, length(fe))
    counts = Vector{Vector{Int}}(undef, length(fe))
    for (d, c) in enumerate(fe)
        remap = zeros(Int32, maximum(c; init = Int32(0)))
        cd = Vector{Int32}(undef, length(rows))
        cnt = Int[]
        @inbounds for (i, r) in enumerate(rows)
            g = c[r]
            k = remap[g]
            if k == 0
                push!(cnt, 0)
                k = remap[g] = Int32(length(cnt))
            end
            cd[i] = k
            cnt[k] += 1
        end
        codes[d] = cd
        counts[d] = cnt
    end
    return _FESample(codes, counts)
end

# v[r] -= mean of v over the group of row r, for one dimension.
function _demean_dim!(
        v::AbstractVector{Float64}, codes::Vector{Int32}, counts::Vector{Int},
        sums::Vector{Float64})
    L = length(counts)
    fill!(view(sums, 1:L), 0.0)
    @inbounds for r in eachindex(v)
        sums[codes[r]] += v[r]
    end
    @inbounds for g in 1:L
        sums[g] /= counts[g]
    end
    @inbounds for r in eachindex(v)
        v[r] -= sums[codes[r]]
    end
    return v
end

function _sweep!(v, fs::_FESample, sums)
    for d in eachindex(fs.codes)
        _demean_dim!(v, fs.codes[d], fs.counts[d], sums)
    end
    return v
end

function _maxabsdiff(a, b)
    m = 0.0
    @inbounds for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end
    return m
end

# Within-transform v in place; returns (sweeps, converged).
function _absorb!(v::AbstractVector{Float64}, fs::_FESample, ::OneSweep, bufs)
    _sweep!(v, fs, bufs.sums)
    return 1, true
end

function _absorb!(v::AbstractVector{Float64}, fs::_FESample, m::AlternatingProjections, bufs)
    if length(fs.codes) == 1
        _demean_dim!(v, fs.codes[1], fs.counts[1], bufs.sums)
        return 1, true
    end
    gx, ggx, sums = bufs.gx, bufs.ggx, bufs.sums
    # The threshold is relative to the column after one sweep, which has removed its level
    # (every dimension spans the constant): the rule depends neither on the units nor on
    # the level of the data, so a small-scale column is not stopped early and an absorbed
    # one still demeans to round-off.
    _sweep!(v, fs, sums)
    thr = m.tol * maximum(abs, v; init = 0.0)
    thr == 0 && return 1, true
    it = 1
    while it < m.maxiter
        copyto!(gx, v)
        _sweep!(gx, fs, sums)
        it += 1
        if _maxabsdiff(gx, v) <= thr
            copyto!(v, gx)
            return it, true
        end
        if !m.accelerate
            copyto!(v, gx)
            continue
        end
        copyto!(ggx, gx)
        _sweep!(ggx, fs, sums)
        it += 1
        if _maxabsdiff(ggx, gx) <= thr
            copyto!(v, ggx)
            return it, true
        end
        # Irons–Tuck: v ← GGv − (Δ₁·Δ₂ / Δ₂·Δ₂) Δ₁, with Δ₁ = GGv − Gv, Δ₂ = GGv − 2Gv + v.
        # Δ₁ lies in the span of the fixed effects, so v stays in v₀ + span and the limit is
        # still the within transform of v₀. Skip the step when Δ₂ is negligible relative
        # to Δ₁ (the ratio, not an absolute size, keeps the step at every scale).
        num = 0.0
        den = 0.0
        n1 = 0.0
        @inbounds for i in eachindex(v)
            d1 = ggx[i] - gx[i]
            d2 = ggx[i] - 2gx[i] + v[i]
            num += d1 * d2
            den += d2 * d2
            n1 += d1 * d1
        end
        if den > eps() * n1
            c = num / den
            @inbounds for i in eachindex(v)
                v[i] = ggx[i] - c * (ggx[i] - gx[i])
            end
        else
            copyto!(v, ggx)
        end
    end
    return it, false
end

function _absorb_buffers(n::Int, fs::_FESample)
    (gx = Vector{Float64}(undef, n),
        ggx = Vector{Float64}(undef, n),
        sums = Vector{Float64}(undef, maximum(length, fs.counts; init = 0)))
end

"""
    _absorb_matrix!(A, fs, method) -> (sweeps, converged)

Within-transform every column of `A` in place. Non-convergence is reported by `panellp`,
once per call, from the fits' `converged` flags.
"""
function _absorb_matrix!(A::AbstractMatrix{Float64}, fs::_FESample, method::AbsorptionMethod)
    isempty(fs.codes) && return 0, true
    bufs = _absorb_buffers(size(A, 1), fs)
    maxit = 0
    conv = true
    for j in axes(A, 2)
        it, c = _absorb!(view(A, :, j), fs, method, bufs)
        maxit = max(maxit, it)
        conv &= c
    end
    return maxit, conv
end

# Connected components of the bipartite graph linking the groups of two dimensions.
function _components(a::Vector{Int32}, b::Vector{Int32}, La::Int, Lb::Int)
    parent = collect(1:(La + Lb))
    function root(x)
        @inbounds while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    @inbounds for r in eachindex(a)
        ra = root(Int(a[r]))
        rb = root(La + Int(b[r]))
        ra != rb && (parent[ra] = rb)
    end
    return count(i -> root(i) == i, 1:(La + Lb))
end

"""
    _fe_components(fs) -> Matrix{Int}

Connected components `m[i, j]` (i < j) of the group graph of every pair of fixed-effect
dimensions on the sample.
"""
function _fe_components(fs::_FESample)
    D = length(fs.codes)
    m = zeros(Int, D, D)
    L = length.(fs.counts)
    for i in 1:D, j in (i + 1):D

        m[i, j] = _components(fs.codes[i], fs.codes[j], L[i], L[j])
    end
    return m
end

"""
    _df_fe(L, m, dims, charge_first) -> Int

Absorbed degrees of freedom of the fixed-effect dimensions `dims` (in formula order), by
reghdfe's rule (reghdfe.mata, `estimate_dof`): every dimension but the first loses
`max_{i<j} m[i, j]` redundant levels (its connected components with each earlier
dimension); the first loses one level when `charge_first` (reghdfe charges it when some
other dimension is nested in a cluster variable). Exact for one or two dimensions; for
three or more it is an upper bound, the conservative count reghdfe reports.
"""
function _df_fe(L::Vector{Int}, m::Matrix{Int}, dims::Vector{Int}, charge_first::Bool)
    isempty(dims) && return 0
    df = L[dims[1]] - (charge_first ? 1 : 0)
    for (k, j) in enumerate(dims)
        k == 1 && continue
        df += L[j] - maximum(m[i, j] for i in dims[1:(k - 1)])
    end
    return df
end

"Whether every group of `fcodes` falls within a single group of `gcodes`."
function _isnested(fcodes::Vector{Int32}, L::Int, gcodes::AbstractVector{Int32})
    seen = zeros(Int32, L)
    @inbounds for r in eachindex(fcodes)
        f = fcodes[r]
        g0 = seen[f]
        if g0 == 0
            seen[f] = gcodes[r]
        elseif g0 != gcodes[r]
            return false
        end
    end
    return true
end
