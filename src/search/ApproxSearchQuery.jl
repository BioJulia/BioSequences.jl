# Approxiamte Search
# ==================
#
# Approximate sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    ApproximateSearchQuery{F<:Function,S<:BioSequence}

Query type for approximate sequence search.

These queries are used as a predicate for the `Base.findnext`, `Base.findprev`,
`Base.occursin`, `Base.findfirst`, and `Base.findlast` functions.

Using these functions with these queries allows you to search a given sequence
for a sub-sequence, whilst allowing a specific number of errors.

In other words they find a subsequence of the target sequence within a specific
[Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of the
query sequence.

# Examples

```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> query = ApproximateSearchQuery(dna"AGGG");

julia> findfirst(query, 0, seq)  # nothing matches with no errors
nothing

julia> findfirst(query, 1, seq)  # seq[3:6] matches with one error
3:6

julia> findfirst(query, 2, seq)  # seq[1:4] matches with two errors
1:4

```

You can pass a comparator function such as `isequal` or `iscompatible` to its
constructor to modify the search behaviour.

The default is `isequal`, however, in biology, sometimes we want a more flexible
comparison to find subsequences of _compatible_ symbols.

```jldoctest
julia> query = ApproximateSearchQuery(dna"AGGG", iscompatible);

julia> occursin(query, 1, dna"AAGNGG")    # 1 mismatch permitted (A vs G) & matched N
true

julia> findnext(query, 1, dna"AAGNGG", 1) # 1 mismatch permitted (A vs G) & matched N
1:4

```

!!! note
    This method of searching for motifs was implemented with smaller query motifs
    in mind.
    
    If you are looking to search for imperfect matches of longer sequences in this
    manner, you are likely better off using some kind of local-alignment algorithm
    or one of the BLAST variants.
"""
struct ApproximateSearchQuery{F<:Function,S<:BioSequence}
    comparator::F           # comparator function
    seq::S                  # query sequence
    fPcom::Vector{UInt64}   # compatibility vector for forward search
    bPcom::Vector{UInt64}   # compatibility vector for backward search
    H::Vector{Int}          # distance vector for alignback function
end

"""
    ApproximateSearchQuery(pat::S, comparator::F = isequal) where {F<:Function,S<:BioSequence}

Construct an [`ApproximateSearchQuery`](@ref) predicate for use with Base find functions.

# Arguments
- `pat`: A concrete BioSequence that is the sub-sequence you want to search for.
- `comparator`: A function used to compare the symbols between sequences. `isequal` by default.
"""
function ApproximateSearchQuery(pat::S, comparator::F = isequal) where {F<:Function,S<:BioSequence}
    m = length(pat)
    if m > 64
        throw(ArgumentError("query pattern sequence must have length of 64 or less"))
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
up to `k` errors.

Symbol comparison is done using the predicate supplied to the query.
By default, `ApproximateSearchQuery`'s predicate is `isequal`.
"""
function Base.findnext(query::ApproximateSearchQuery, k::Integer, seq::BioSequence, start::Integer)
    return _approxsearch(query, k, seq, start, lastindex(seq), true)
end

"""
    findprev(query, k, seq, start)

Return the range of the last occurrence of `query` in `seq[stop:start]` allowing
up to `k` errors.

Symbol comparison is done using the predicate supplied to the query.
By default, `ApproximateSearchQuery`'s predicate is `isequal`.
"""
function Base.findprev(query::ApproximateSearchQuery, k::Integer, seq::BioSequence, start::Integer)
    return _approxsearch(query, k, seq, start, firstindex(seq), false)
end

Base.findfirst(query::ApproximateSearchQuery, k::Integer, seq::BioSequence) = findnext(query, k, seq, firstindex(seq))
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
