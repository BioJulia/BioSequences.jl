"""
    ExactSearchQuery{F<:Function,S<:BioSequence}

Query type for exact sequence search.

An exact search, is one where are you are looking in some given sequence, for
exact instances of some given substring.

These queries are used as a predicate for the `Base.findnext`, `Base.findprev`,
`Base.occursin`, `Base.findfirst`, and `Base.findlast` functions.

# Examples

```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> query = ExactSearchQuery(dna"AGC");

julia> findfirst(query, seq)
3:5

julia> findlast(query, seq)
8:10

julia> occursin(query, seq)
true

```

You can pass a comparator function such as `isequal` or `iscompatible` to its
constructor to modify the search behaviour.

The default is `isequal`, however, in biology, sometimes we want a more flexible
comparison to find subsequences of _compatible_ symbols.

```jldoctest
julia> findfirst(ExactSearchQuery(dna"CGT", iscompatible), dna"ACNT")  # 'N' matches 'G'
2:4

julia> findfirst(ExactSearchQuery(dna"CNT", iscompatible), dna"ACGT")  # 'G' matches 'N'
2:4

julia> occursin(ExactSearchQuery(dna"CNT", iscompatible), dna"ACNT")
true

```
"""
struct ExactSearchQuery{F<:Function,S<:BioSequence}
    comparator::F       # comparator function
    seq::S              # query sequence
    bloom_mask::UInt64  # compatibility bits / bloom mask
    fshift::Int         # shift length for forward search
    bshift::Int         # shift length for backward search
end

"""
    ExactSearchQuery(pat::BioSequence, comparator::Function = isequal)

Construct an [`ExactSearchQuery`](@ref) predicate for use with Base find functions.

# Arguments
- `pat`: A concrete BioSequence that is the sub-sequence you want to search for.
- `comparator`: A function used to compare the symbols between sequences. `isequal` by default.
"""
function ExactSearchQuery(pat::BioSequence, comparator::Function = isequal)
    T = ExactSearchQuery{typeof(comparator),typeof(pat)}
    
    m = length(pat)
    if m == 0
        return T(comparator, pat, UInt64(0), 0, 0)
    end
    
    first = pat[1]
    last = pat[end]
    bloom_mask = zero(UInt64)
    fshift = bshift = m
    
    for i in 1:lastindex(pat)
        x = pat[i]
        bloom_mask |= _bloom_bits(typeof(comparator), x)
        if comparator(x, last) && i < m
            fshift = m - i
        end
    end
    
    for i in lastindex(pat):-1:1
        x = pat[i]
        if comparator(x, first) && i > 1
            bshift = i - 1
        end
    end
    
    return T(comparator, pat, bloom_mask, fshift, bshift)
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

Symbol comparison is done using the predicate supplied to the query.
By default, `ExactSearchQuery`'s predicate is `isequal`.
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

