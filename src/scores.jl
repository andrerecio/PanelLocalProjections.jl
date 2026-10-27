# Sufficient statistics of the scores of one horizon. For a variance computed after
# estimation we need, for each clustering dimension, the score sums by group; for two-way
# clustering also the sums by intersection group; and Σ u u' for the heteroskedasticity
# terms. Group codes are global (the same group has the same code at every horizon), so
# cross-horizon covariances can be formed from the stored sums as well.

# Score sums by group: S[g, :] = Σ_{r: codes[r] = g} U[r, :] over the kept columns.
function _groupsums(U::AbstractMatrix{Float64}, codes::AbstractVector{Int32}, G::Int,
        kidx::Vector{Int}, K::Int)
    S = zeros(G, K)
    @inbounds for (c, j) in enumerate(kidx)
        for r in axes(U, 1)
            S[codes[r], j] += U[r, c]
        end
    end
    return S
end

_npresent(codes::AbstractVector{Int32}, G::Int) = count(_present(codes, G))
function _present(codes::AbstractVector{Int32}, G::Int)
    seen = falses(G)
    @inbounds for c in codes
        seen[c] = true
    end
    return seen
end

# How the intersection of dimensions a and b relates to them on the sample:
# :obs (every intersection group is one row), :first (a nested in b, so a ∩ b = a),
# :second (b nested in a), or :own (a genuine coarser intersection, stored).
function _pairkind(ca::AbstractVector{Int32}, Ga::Int, cb::AbstractVector{Int32}, Gb::Int)
    _isnested(Vector{Int32}(ca), Ga, cb) && return :first
    _isnested(Vector{Int32}(cb), Gb, ca) && return :second
    d = Dict{Tuple{Int32, Int32}, Int}()
    @inbounds for r in eachindex(ca)
        k = (ca[r], cb[r])
        d[k] = get(d, k, 0) + 1
    end
    return all(==(1), values(d)) ? :obs : :own
end

function _pairsums(U, ca, cb, kidx, K)
    d = Dict{Tuple{Int32, Int32}, Int32}()
    codes = Vector{Int32}(undef, length(ca))
    @inbounds for r in eachindex(ca)
        k = (ca[r], cb[r])
        g = get(d, k, Int32(0))
        g == 0 && (g = d[k] = Int32(length(d) + 1))
        codes[r] = g
    end
    return _groupsums(U, codes, length(d), kidx, K)
end

"""
    _score_stats(U, keep, groups, ngroups) -> (ScoreStats, groups present)

`U` holds the scores u_r = d_r ê_r on the kept columns (rows = estimation sample);
`groups[g]` are the codes of score dimension g on the same rows.
"""
function _score_stats(U::AbstractMatrix{Float64}, keep::BitVector,
        groups::Vector{<:AbstractVector{Int32}}, ngroups::Vector{Int})
    K = length(keep)
    kidx = findall(keep)
    Mobs = zeros(K, K)
    Mobs[kidx, kidx] = U' * U
    S = [_groupsums(U, groups[g], ngroups[g], kidx, K) for g in eachindex(groups)]
    present = [_npresent(groups[g], ngroups[g]) for g in eachindex(groups)]
    pairs = Tuple{Int, Int}[]
    Spair = Matrix{Float64}[]
    kinds = Symbol[]
    for a in eachindex(groups), b in (a + 1):length(groups)

        push!(pairs, (a, b))
        # unit × time is the observation level by construction
        kind = (a, b) == (1, 2) ? :obs :
               _pairkind(groups[a], ngroups[a], groups[b], ngroups[b])
        push!(kinds, kind)
        push!(Spair, kind === :own ? _pairsums(U, groups[a], groups[b], kidx, K) :
                     zeros(0, K))
    end
    return ScoreStats(Mobs, S, pairs, Spair, kinds), present
end
