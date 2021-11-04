
# Exact Search
# ============
#
# Exact sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

abstract type ExactSearchQuery{S<:BioSequence} end

@inline function checkeltype(seq1::BioSequence, seq2::BioSequence)
    if eltype(seq1) != eltype(seq2)
        throw(ArgumentError("the element type of two sequences must match"))
    end
end

@inline function _check_bloom_mask(q::ExactSearchQuery, x::BioSymbol)
    return _bloom_mask(q) & _bloom_bits(typeof(q), x) == 0
end


"""
    findnext(query::ExactSearchQuery, seq::BioSequence, start::Integer)

Return the index of the first occurrence of `query` in `seq`.

Symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findnext(query::ExactSearchQuery, seq::BioSequence, start::Integer)
    checkeltype(seq, query.seq)
    i = quicksearch(query, seq, start, lastindex(seq))
    if i == 0
        return nothing
    else
        return i:i+length(query.seq)-1
    end
end

Base.findfirst(pat::ExactSearchQuery, seq::BioSequence) = findnext(pat, seq, firstindex(seq))


"""
    findprev(query::ExactSearchQuery, seq::BioSequence, start::Integer)

Return the index of the last occurrence of `query` in `seq`.
Symbol comparison is done using `BioSequences.iscompatible`.
"""
function Base.findprev(query::ExactSearchQuery, seq::BioSequence, start::Integer)
    checkeltype(seq, query.seq)
    i = quickrsearch(query, seq, start, 1)
    if i == 0
        return nothing
    else
        return i:i+length(query.seq)-1
    end
end

Base.findlast(query::ExactSearchQuery, seq::BioSequence) = findprev(query, seq, lastindex(seq))


"""
    occursin(x::MatchSearchQuery, y::BioSequence)

Return Bool indicating presence of exact match of x in y.
"""
Base.occursin(x::MatchSearchQuery, y::BioSequence) = quicksearch(x, y, 1, lastindex(y)) != 0


"""
    query_preprocess(type, query)

Preprocesses a search query by building the bloom mask
and computing shift values for quicksearch algorithm in
advance.
"""
function query_preprocess(type::Type{<:ExactSearchQuery{S}}, query::S) where {S<:BioSequence}
    if length(query) == 0
        return UInt64(0), 0, 0
    end
    m = length(query)
    first = query[1]
    last = query[end]
    bloom_mask = zero(UInt64)
    fshift = bshift = m
    for i in 1:lastindex(query)
        x = query[i]
        bloom_mask |= _bloom_bits(type, x)
        if _ismatch(type, x, last) && i < m
            fshift = m - i
        end
    end
    for i in lastindex(query):-1:1
        x = query[i]
        if _ismatch(type, x, first) && i > 1
            bshift = i - 1
        end
    end
    return bloom_mask, fshift, bshift
end


"""
Query type for exact sequence search.
"""
struct MatchSearchQuery{S<:BioSequence} <: ExactSearchQuery{S}
    seq::S         # query sequence
    bloom_mask::UInt64  # compatibility bits / bloom mask
    fshift::Int    # shift length for forward search
    bshift::Int    # shift length for backward search
end

function MatchSearchQuery(query::BioSequence)
    bf, fs, bs = query_preprocess(MatchSearchQuery{typeof(query)}, query)
    return MatchSearchQuery(query, bf, fs, bs)
end

@inline _bloom_mask(q::MatchSearchQuery) = q.bloom_mask

@inline function _bloom_bits(::Type{<:MatchSearchQuery}, x::BioSymbol)
    return (UInt64(1) << (encoded_data(x) & 63))
end

@inline _ismatch(::Type{<:MatchSearchQuery}, x::BioSymbol, y::BioSymbol) = x == y

# This algorithm is borrowed from base/strings/search.jl, which looks similar to
# Sunday's Quick Search algorithm.
function quicksearch(query::MatchSearchQuery, seq, start, stop)
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
        if pat[m] == seq[s + m]
            # If last element in pattern matches...
            i = m - 1
            while i > 0 # Consider previous element in pattern...
                if pat[i] != seq[s + i]
                    # If it does not match, break.
                    break
                end
                i -= 1 # Else, consider next previous element in pattern...
            end
            if i == 0 # If I is 0, every element in pattern matched.
                return s + 1  # found
            # No match, attempt to rule out next element if end not reached...
            # Essentially, if the next base 
            elseif s < stop′ && _check_bloom_mask(query, seq[s + m + 1])
                s += m + 1
            else
                s += query.fshift
            end
        elseif s < stop′ && _check_bloom_mask(query, seq[s + m + 1])
            s += m + 1
        else
            s += 1
        end
    end

    return 0  # not found
end


function quickrsearch(query::MatchSearchQuery, seq, start, stop)
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
        if pat[1] == seq[s + 1]
            i = 2
            while i < m + 1
                if pat[i] != seq[s + i]
                    break
                end
                i += 1
            end
            if i == m + 1
                return s + 1  # found
            elseif s > stop′ && _check_bloom_mask(query, seq[s])
                s -= m + 1
            else
                s -= query.bshift
            end
        elseif s > stop′ && _check_bloom_mask(query, seq[s])
            s -= m + 1
        else
            s -= 1
        end
    end

    return 0  # not found
end


"""
Query type for compatible sequence search.
"""
struct CompatibleSearchQuery{S<:BioSequence}
    seq::S         # query sequence
    cbits::UInt32  # compatibility bits
    fshift::Int    # shift length for forward search
    bshift::Int    # shift length for backward search
end


function CompatibleSearchQuery(query::BioSequence)
    cbits, fshift, bshift = search_preprocess(query)
    return CompatibleSearchQuery(query, cbits, fshift, bshift)
end

@inline bloom_mask(q::CompatibleSearchQuery) = q.bloom_mask
@inline bloom_bits(::Type{<:CompatibleSearchQuery}, x::NucleicAcid) = encoded_data(x)

@inline function ismatch(::Type{<:CompatibleSearchQuery}, x::BioSymbol, y::BioSymbol)
    return iscompatible(x, y)
end

function quicksearch(query::CompatibleSearchQuery, seq, start, stop)
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
        if iscompatible(pat[m], seq[s+m])
            i = m - 1
            while i > 0
                if !iscompatible(pat[i], seq[s+i])
                    break
                end
                i -= 1
            end
            if i == 0
                return s + 1  # found
            elseif s < stop′ && (query.cbits & compatbits(seq[s+m+1]) == 0)
                s += m + 1
            elseif isambiguous(seq[s+m])
                s += 1
            else
                s += query.fshift
            end
        elseif s < stop′ && (query.cbits & compatbits(seq[s+m+1]) == 0)
            s += m + 1
        else
            s += 1
        end
    end

    return 0  # not found
end


function quickrsearch(query::CompatibleSearchQuery, seq, start, stop)
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
        if iscompatible(pat[1], seq[s+1])
            i = 2
            while i < m + 1
                if !iscompatible(pat[i], seq[s+i])
                    break
                end
                i += 1
            end
            if i == m + 1
                return s + 1  # found
            elseif s > stop′ && (query.cbits & compatbits(seq[s]) == 0)
                s -= m + 1
            elseif isambiguous(seq[s+1])
                s -= 1
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
