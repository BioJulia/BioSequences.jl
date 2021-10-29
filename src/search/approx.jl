# Approxiamte Search
# ==================
#
# Approximate sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Query type for approximate sequence search.
"""
struct ApproximateSearchQuery{S<:BioSequence}
    seq::S          # query sequence
    k::Int          # number of mismatches permitted
    fPcom::Vector   # compatibility vector for forward search
    bPcom::Vector   # compatibility vector for backward search
    H::Vector{Int}  # distance vector for alignback function

    function ApproximateSearchQuery{S}(seq::BioSequence, k::Integer, direction::Symbol) where S
        if direction == :forward
            fPcom = approx_preprocess(seq, true)
            bPcom = []
        elseif direction == :backward
            fPcom = []
            bPcom = approx_preprocess(seq, false)
        elseif direction == :both
            fPcom = approx_preprocess(seq, true)
            bPcom = approx_preprocess(seq, false)
        else
            throw(ArgumentError("direction '$direction' is invalid"))
        end
        H = Vector{Int}(undef, length(seq) + 1)
        return new{S}(seq, k, fPcom, bPcom, H)
    end
end

Base.isempty(x::ApproximateSearchQuery) = isempty(x.seq)

"""
    ApproximateSearchQuery(pat::BioSequence, k::Integer[, direction=:both])

Create an query object for approximate sequence search from the `pat` sequence.

# Arguments
* `pat`: Query sequence.
* `k`: The number of mismatches /errors permitted.
* `direction=:both`: Search direction (`:forward`, `:backward`, or `:both`).
"""
function ApproximateSearchQuery(pat::BioSequence, k::Integer, direction::Symbol = :both)
    return ApproximateSearchQuery{typeof(pat)}(pat, k, direction)
end

function approx_preprocess(pat, forward)
    # select a bit vector type
    # TODO: BigInt is very slow, consider implementing "4.2 THE BLOCKS MODEL"
    m = length(pat)
    T = m ≤ 64 ? UInt64 : m ≤ 128 ? UInt128 : BigInt
    Σ = alphabet(eltype(pat))
    Pcom = zeros(T, length(Σ))
    for i in 1:m
        y = forward ? pat[i] : pat[end-i+1]
        for x in Σ
            if BioSequences.iscompatible(x, y)
                Pcom[reinterpret(UInt8, x)+0x01] |= one(T) << (i - 1)
            end
        end
    end
    return Pcom
end

"""
    findnext(query, seq[, start=1])

Return the range of the first occurrence of `pat` in `seq[start:stop]` allowing
up to `k` errors; symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findnext(query::ApproximateSearchQuery, seq::BioSequence,
                       start::Integer = firstindex(seq))
    return _approxsearch(query, seq, k, start, lastindex(seq), true)
end

# TODO: Needed?
Base.findfirst(query::ApproximateSearchQuery, seq::BioSequence) = findnext(query, seq)

"""
    findprev(query, seq, k[, start=lastindex(seq)[, stop=1]])

Return the range of the last occurrence of `pat` in `seq[stop:start]` allowing
up to `k` errors; symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findprev(query::ApproximateSearchQuery, seq::BioSequence,
                       start::Integer = lastindex(seq))
    return _approxsearch(query, seq, k, start, firstindex(seq), false)
end

# TODO: Needed?
Base.findlast(query::ApproximateSearchQuery, seq::BioSequence) = findprev(query, seq)

function _approxsearch(query, seq, start, stop, forward)
    if forward && isempty(query.fPcom)
        throw(ArgumentError("query is not preprocessed for forward search"))
    end
    if !forward && isempty(query.bPcom)
        throw(ArgumentError("query is not preprocessed for backward search"))
    end

    if k ≥ length(query.seq)
        return nothing
    end

    # search the approximate suffix
    sas = search_approx_suffix(
        forward ? query.fPcom : query.bPcom,
        query.seq, seq, query.k, start, stop, forward)
    if matchstop == 0
        return 0:-1
    end
    matchstop, dist = sas

    # locate the starting position of the match
    matchstart = alignback!(query.H, query.seq, seq, dist, start, matchstop, forward)
    if forward
        return matchstart:matchstop
    else
        return matchstop:matchstart
    end
end

# This returns the end index of a suffix sequence with up to `k` errors.
# More formally, when `forward = true`, it returns the minimum `j ∈ start:stop`
# such that `min_{g ∈ 1:j} δ(pat, seq[g:j]) ≤ k` where `δ(s, t)` is the edit
# distance between `s` and `t` sequences. See Myers' paper for details:
# Myers, Gene. "A fast bit-vector algorithm for approximate string matching
# based on dynamic programming." Journal of the ACM (JACM) 46.3 (1999): 395-415.
# NOTE: `Pcom` corresponds to `Peq` in the paper.
function search_approx_suffix(Pcom::Vector{T}, pat, seq, k, start, stop, forward) where {T}
    if k < 0
        throw(ArgumentError("the number of errors must be non-negative"))
    end

    m = length(pat)
    n = length(seq)
    @assert T == BigInt || m ≤ sizeof(T) * 8

    Pv::T = (one(T) << m) - one(T)
    Mv::T = zero(T)
    dist = m
    j = forward ? max(start, 1) : min(start, n)

    if dist ≤ k
        return j, dist
    end

    while (forward && j ≤ min(stop, n)) || (!forward && j ≥ max(stop, 1))
        Eq = Pcom[reinterpret(UInt8, seq[j])+0x01]
        Xv = Eq | Mv
        Xh = (((Eq & Pv) + Pv) ⊻ Pv) | Eq

        Ph = Mv | ~(Xh | Pv)
        Mh = Pv & Xh
        if (Ph >> (m - 1)) & 1 != 0
            dist += 1
        elseif (Mh >> (m - 1)) & 1 != 0
            dist -= 1
        end

        if dist ≤ k
            return j, dist  # found
        end

        Ph <<= 1
        Mh <<= 1
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv
        j += ifelse(forward, +1, -1)
    end

    return nothing  # not found
end

# run dynamic programming to get the starting position of the alignment
function alignback!(H, pat, seq, dist, start, matchstop, forward)
    m = length(pat)
    n = length(seq)

    # initialize the cost column
    for i in 0:m
        H[i+1] = i
    end

    j = ret = matchstop
    found = false
    while (forward && j ≥ max(start, 1)) || (!forward && j ≤ min(start, n))
        y = seq[j]
        h_diag = H[1]
        for i in 1:m
            x = forward ? pat[end-i+1] : pat[i]
            h = min(
                H[i] + 1,
                H[i+1] + 1,
                h_diag + ifelse(iscompatible(x, y), 0, 1))
            h_diag = H[i+1]
            H[i+1] = h
        end
        if H[m+1] == dist
            ret = j
            found = true
        end
        j += ifelse(forward, -1, +1)
    end
    @assert found

    return ret
end
