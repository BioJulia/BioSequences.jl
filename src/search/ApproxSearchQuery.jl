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
struct ApproximateSearchQuery{F<:Function,S<:BioSequence}
    comparator::F           # comparator function
    seq::S                  # query sequence
    fPcom::Vector{UInt64}   # compatibility vector for forward search
    bPcom::Vector{UInt64}   # compatibility vector for backward search
    H::Vector{Int}          # distance vector for alignback function
end

function ApproximateSearchQuery(pat::S, comparator::F = iscompatible) where {F<:Function,S<:BioSequence}
    m = length(pat)
    if m > 64
        throw(ArgumentEror("query pattern sequence must have length of 64 or less"))
    end
    Σ = alphabet(eltype(pat))
    fw = zeros(UInt64, length(Σ))
    for i in 1:m
        y = pat[i]
        for x in Σ
            if comparator(x, y)
                fw[encoded_data(x) + 0x01] |= UInt(1) << (i - 1)
            end
        end
    end
    shift = 64 - m
    bw = [bitreverse(i) >>> (shift & 63) for i in fw]
    
    H = Vector{Int}(undef, length(pat) + 1)
    
    return ApproximateSearchQuery{F,S}(comparator, copy(pat), fw, bw, H)
end

"""
    findnext(query, k, seq, start)

Return the range of the first occurrence of `pat` in `seq[start:stop]` allowing
up to `k` errors; symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findnext(query::ApproximateSearchQuery, k::Integer, seq::BioSequence, start::Integer)
    return _approxsearch(query, k, seq, start, lastindex(seq), true)
end

# TODO: Needed?
Base.findfirst(query::ApproximateSearchQuery, k::Integer, seq::BioSequence) = findnext(query, k, seq, firstindex(seq))

"""
    findprev(query, k, seq, start)

Return the range of the last occurrence of `query` in `seq[stop:start]` allowing
up to `k` errors; symbol comparison is done using the predicate function supplied to the query. 
By default, `ApproxSearchQuery`'s predicate is`BioSequences.iscompatible`.
"""
function Base.findprev(query::ApproximateSearchQuery, k::Integer, seq::BioSequence, start::Integer)
    return _approxsearch(query, k, seq, start, firstindex(seq), false)
end

# TODO: Needed?
Base.findlast(query::ApproximateSearchQuery, k::Integer, seq::BioSequence) = findprev(query, k, seq, lastindex(seq))

Base.occursin(query::ApproximateSearchQuery, k::Integer, seq::BioSequence) = !isnothing(_approxsearch(query, k, seq, 1, lastindex(seq), true))

function _approxsearch(query::ApproximateSearchQuery, k::Integer, seq::BioSequence, start::Integer, stop::Integer, forward::Bool)
    checkeltype(query.seq, seq)
    if k ≥ length(query.seq)
        return start:start-1
    end

    # search the approximate suffix
    matchstop, dist = search_approx_suffix(
        forward ? query.fPcom : query.bPcom,
        query.seq, seq, k, start, stop, forward)
    if matchstop == 0 # No match
        #return 0:-1
        return nothing
    end

    # locate the starting position of the match
    matchstart = alignback!(query, seq, dist, start, matchstop, forward)
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
function search_approx_suffix(Pcom::Vector{UInt64}, pat::BioSequence, seq::BioSequence, k::Integer, start::Integer, stop::Integer, forward::Bool)
    if k < 0
        throw(ArgumentError("the number of errors must be non-negative"))
    end
    
    m = length(pat)
    n = length(seq)

    Pv::UInt64 = (one(UInt64) << m) - one(UInt64)
    Mv::UInt64 = zero(UInt64)
    dist = m
    j = forward ? max(start, 1) : min(start, n)

    if dist ≤ k
        return j, dist
    end

    while (forward && j ≤ min(stop, n)) || (!forward && j ≥ max(stop, 1))
        Eq = Pcom[reinterpret(UInt8, seq[j]) + 0x01]
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

    return 0, -1  # not found
end

# run dynamic programming to get the starting position of the alignment
function alignback!(query::ApproximateSearchQuery{<:Function,<:BioSequence}, seq::BioSequence, dist::Int, start::Integer, matchstop::Integer, forward::Bool)
    comparator = query.comparator
    H = query.H
    pat = query.seq
    
    m = length(pat)
    n = length(seq)
    
    # initialize the cost column
    for i in 0:m
        H[i + 1] = i
    end
    
    j = ret = matchstop
    found = false
    while (forward && j ≥ max(start, 1)) || (!forward && j ≤ min(start, n))
        y = seq[j]
        h_diag = H[1]
        for i in 1:m
            x = forward ? pat[end - i + 1] : pat[i]
            h = min(
                H[i] + 1,
                H[i + 1] + 1,
                h_diag + ifelse(comparator(x, y), 0, 1))
            h_diag = H[i + 1]
            H[i + 1] = h
        end
        if H[m + 1] == dist
            ret = j
            found = true
        end
        j += ifelse(forward, -1, +1)
    end
    @assert found

    return ret
end
