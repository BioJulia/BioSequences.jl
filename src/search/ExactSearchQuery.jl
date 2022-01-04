"""
Query type for exact sequence search.
"""
struct ExactSearchQuery{F<:Function,S<:BioSequence}
    comparator::F       # comparator function
    seq::S              # query sequence
    bloom_mask::UInt64  # compatibility bits / bloom mask
    fshift::Int         # shift length for forward search
    bshift::Int         # shift length for backward search
end

function ExactSearchQuery(query::BioSequence, comparator::Function = isequal)
    T = ExactSearchQuery{typeof(comparator),typeof(query)}
    
    m = length(query)
    if m == 0
        return T(comparator, query, UInt64(0), 0, 0)
    end
    
    first = query[1]
    last = query[end]
    bloom_mask = zero(UInt64)
    fshift = bshift = m
    
    for i in 1:lastindex(query)
        x = query[i]
        bloom_mask |= _bloom_bits(typeof(comparator), x)
        if comparator(x, last) && i < m
            fshift = m - i
        end
    end
    
    for i in lastindex(query):-1:1
        x = query[i]
        if comparator(x, first) && i > 1
            bshift = i - 1
        end
    end
    
    return T(comparator, query, bloom_mask, fshift, bshift)
end


@inline _check_ambiguous(q::ExactSearchQuery{typeof(isequal),<:BioSequence}) = false
@inline _check_ambiguous(q::ExactSearchQuery{typeof(iscompatible),<:BioSequence}) = true

@inline function _bloom_bits(::Type{typeof(isequal)}, x::BioSymbol)
    return (UInt64(1) << (encoded_data(x) & 63))
end
@inline function _bloom_bits(::Type{typeof(iscompatible)}, x::BioSymbol)
    return compatbits(x)
end

@inline function checkeltype(seq1::BioSequence, seq2::BioSequence)
    if eltype(seq1) != eltype(seq2)
        throw(ArgumentError("the element type of two sequences must match"))
    end
end

function quicksearch(query::ExactSearchQuery, seq::BioSequence, start::Integer, stop::Integer)
    pat = query.seq
    comparator = query.comparator
    bloom_mask = query.bloom_mask
    ambig_check = _check_ambiguous(query)
    
    checkeltype(seq, pat)
    
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
        if comparator(pat[m], seq[s + m])
            i = m - 1
            while i > 0
                if !comparator(pat[i], seq[s + i])
                    break
                end
                i -= 1
            end
            if i == 0
                return s + 1  # found
            elseif s < stop′ && (bloom_mask & _bloom_bits(typeof(comparator), seq[s + m + 1]) == 0)
                s += m + 1
            elseif ambig_check && isambiguous(seq[s + m])
                s += 1
            else
                s += query.fshift
            end
        elseif s < stop′ && (bloom_mask & _bloom_bits(typeof(comparator), seq[s + m + 1]) == 0)
            s += m + 1
        else
            s += 1
        end
    end

    return 0  # not found
end

function quickrsearch(query::ExactSearchQuery, seq::BioSequence, start::Integer, stop::Integer)
    pat = query.seq
    comparator = query.comparator
    bloom_mask = query.bloom_mask
    ambig_check = _check_ambiguous(query)
    
    checkeltype(seq, pat)
    
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
        if comparator(pat[1], seq[s + 1])
            i = 2
            while i < m + 1
                if !comparator(pat[i], seq[s + i])
                    break
                end
                i += 1
            end
            if i == m + 1
                return s + 1  # found
            elseif s > stop′ && (bloom_mask & _bloom_bits(typeof(comparator), seq[s]) == 0)
                s -= m + 1
            elseif ambig_check && isambiguous(seq[s + 1])
                s -= 1
            else
                s -= query.bshift
            end
        elseif s > stop′ && (bloom_mask & _bloom_bits(typeof(comparator), seq[s]) == 0)
            s -= m + 1
        else
            s -= 1
        end
    end

    return 0  # not found
end


"""
    findnext(query::ExactSearchQuery, seq::BioSequence, start::Integer)

Return the index of the first occurrence of `query` in `seq`.

Symbol comparison is done using the predicate supplied to the query.
By default, `ExactSearchQuery`'s predicate is `isequal`.
"""
function Base.findnext(query::ExactSearchQuery, seq::BioSequence, start::Integer)
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
    i = quickrsearch(query, seq, start, 1)
    if i == 0
        return nothing
    else
        return i:i+length(query.seq)-1
    end
end

Base.findlast(query::ExactSearchQuery, seq::BioSequence) = findprev(query, seq, lastindex(seq))

"""
    occursin(x::ExactSearchQuery, y::BioSequence)

Return Bool indicating presence of exact match of x in y.
"""
Base.occursin(x::ExactSearchQuery, y::BioSequence) = quicksearch(x, y, 1, lastindex(y)) != 0

