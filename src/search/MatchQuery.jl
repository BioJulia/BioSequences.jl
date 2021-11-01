
# Exact Search
# ============
#
# Exact sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Reference Implementation
# ------------------------
#
# `searchindex` and `rsearchindex` should work exactly the same way as the
# following reference implementation:
#=
# Return the index of the first occurrence of `pat` in `seq[start:stop]`.
function searchindex(seq, pat, start=1, stop=lastindex(seq))
    m = length(pat)
    n = length(seq)
    for s in max(start-1, 0):min(stop, n)-m
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

# Return the index of the last occurrence of `pat` in `seq[stop:start]`.
function rsearchindex(seq, pat, start=lastindex(seq), stop=1)
    n = length(seq)
    m = length(pat)
    for s in min(start, n)-m:-1:max(stop-1, 0)
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

function occurs_with_shift(pat, seq, shift)
    for i in 1:lastindex(pat)
        if !iscompatible(pat[i], seq[shift+i])
            return false
        end
    end
    return true
end
=#


# Actual Implementation
# ---------------------

"""
Query type for exact sequence search.
"""
struct MatchQuery{S<:BioSequence}
    seq::S         # query sequence
    cbits::UInt32  # compatibility bits / bloom mask
    fshift::Int    # shift length for forward search
    bshift::Int    # shift length for backward search
end

function MatchQuery(query::BioSequence)
    cbits, fshift, bshift = preprocess(query)
    return MatchQuery(query, cbits, fshift, bshift)
end

"""
    preprocess(query::MatchQuery)

Preprocesses a search query by building the bloom mask
and computing shift values for quicksearch algorithm in
advance.
"""
function preprocess(query::MatchQuery)
    if length(query) == 0
        return UInt32(0), 0, 0
    end
    m = length(query)
    first = query[1]
    last = query[end]
    cbits::UInt32 = 0
    fshift = bshift = m
    for i in 1:lastindex(query)
        x = query[i]
        cbits |= compatbits(x)
        if x == last && i < m
            fshift = m - i
        end
    end
    for i in lastindex(query):-1:1
        x = query[i]
        if x == first && i > 1
            bshift = i - 1
        end
    end
    return cbits, fshift, bshift
end

function checkeltype(seq1, seq2)
    if eltype(seq1) != eltype(seq2)
        throw(ArgumentError("the element type of two sequences must match"))
    end
end


# Forward
# -------

# This algorithm is borrowed from base/strings/search.jl, which looks similar to
# Sunday's Quick Search algorithm.
function quicksearch(query::MatchQuery, seq, start, stop)
    pat = query.seq
    m = length(pat)
    n = length(seq)
    stop′ = min(stop, n) - m
    s::Int = max(start - 1, 0)

    if m == 0  # empty query
        if s ≤ stop′
            return s + 1  # found
        else
            return 0  # not found
        end
    end

    while s ≤ stop′
        if pat[m] == seq[s + m] # If last element in pattern matches...
            i = m - 1
            while i > 0 # Consider previous element in pattern...
                if pat[i] != seq[s + i] # If it does not match, break.
                    break
                end
                i -= 1 # Else, consider next previous element in pattern...
            end
            if i == 0 # If I is 0, every element in pattern matched.
                return s + 1  # found
            # No match, attempt to rule out next element if end not reached...
            # Essentially, if the next base 
            elseif s < stop′ && (query.cbits & compatbits(seq[s + m + 1]) == 0)
                s += m + 1
            #elseif isambiguous(seq[s+m])
            #    s += 1
            else
                s += query.fshift
            end
        elseif s < stop′ && (query.cbits & compatbits(seq[s + m + 1]) == 0)
            s += m + 1
        else
            s += 1
        end
    end

    return 0  # not found
end

"""
    findnext(query::MatchQuery, seq::BioSequence, start::Integer)

Return the index of the first occurrence of `query` in `seq`.

Symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findnext(query::MatchQuery, seq::BioSequence, start::Integer)
    checkeltype(seq, query.seq)
    i = quicksearch(query, seq, start, lastindex(seq))
    if i == 0
        return nothing
    else
        return i:i+length(query.seq)-1
    end
end

Base.findfirst(pat::MatchQuery, seq::BioSequence) = findnext(pat, seq, firstindex(seq))


# Backward
# --------

function quickrsearch(query::MatchQuery, seq, start, stop)
    pat = query.seq
    m = length(pat)
    n = length(seq)
    stop′ = max(stop - 1, 0)
    s::Int = min(start, n) - m

    if m == 0  # empty query
        if s ≥ stop′
            return s + 1  # found
        else
            return 0  # not found
        end
    end

    while s ≥ stop′
        if pat[1] == seq[s+1]
            i = 2
            while i < m + 1
                if pat[i] != seq[s+i]
                    break
                end
                i += 1
            end
            if i == m + 1
                return s + 1  # found
            elseif s > stop′ && (query.cbits & compatbits(seq[s]) == 0)
                s -= m + 1
            #elseif isambiguous(seq[s+1])
            #    s -= 1
            else
                s -= query.bshift
            end
        elseif s > stop′ && (query.cbits & compatbits(seq[s]) == 0)
            s -= m + 1
        else
            s -= 1
        end
    end

    return 0  # not found
end

"""
    findprev(query::MatchQuery, seq::BioSequence, start::Integer)

Return the index of the last occurrence of `query` in `seq`.
Symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findprev(query::MatchQuery, seq::BioSequence, start::Integer)
    checkeltype(seq, query.seq)
    i = quickrsearch(seq, query, start, 1)
    if i == 0
        return nothing
    else
        return i:i+length(query.seq)-1
    end
end

Base.findlast(query::MatchQuery, seq::BioSequence) = findprev(query, seq, lastindex(seq))

"""
    occursin(x::MatchQuery, y::BioSequence)

Return Bool indicating presence of exact match of x in y.
"""
Base.occursin(x::MatchQuery, y::BioSequence) = quicksearch(x, y, 1, lastindex(y)) != 0
