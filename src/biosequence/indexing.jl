###
### Indexing
###
###
### Indexing methods for mutable biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

@inline function Base.iterate(seq::BioSequence, i::Int = firstindex(seq))
    i > lastindex(seq) ? nothing : @inbounds seq[i], i + 1
end

## Bounds checking
function Base.checkbounds(x::BioSequence, i::Integer)
    firstindex(x) ≤ i ≤ lastindex(x) || throw(BoundsError(x, i))
end

function Base.checkbounds(x::BioSequence, locs::AbstractVector{Bool})
    length(x) == length(locs) || throw(BoundsError(x, lastindex(locs)))
end

@inline function Base.checkbounds(seq::BioSequence, locs::AbstractVector)
    for i in locs
        checkbounds(seq, i)
    end
    return true
end

@inline function Base.checkbounds(seq::BioSequence, range::UnitRange)
    if !isempty(range) && (first(range) < 1 || last(range) > length(seq))
        throw(BoundsError(seq, range))
    end
end

## Getindex
Base.@propagate_inbounds function Base.getindex(x::BioSequence, i::Integer)
    @boundscheck checkbounds(x, i)
    data = extract_encoded_element(x, i)
    return decode(Alphabet(x), data)
end

Base.@propagate_inbounds function Base.getindex(x::BioSequence, bools::AbstractVector{Bool})
    @boundscheck checkbounds(x, bools)
    res = typeof(x)(undef, count(bools))
    ind = 0
    @inbounds for i in eachindex(bools)
        if bools[i]
            ind += 1
            res[ind] = x[i]
        end
    end
    res
end

Base.@propagate_inbounds function Base.getindex(x::BioSequence, i::AbstractVector{<:Integer})
    @boundscheck checkbounds(x, i)
    isempty(i) && return empty(x)
    res = typeof(x)(undef, length(i))
    @inbounds for ind in eachindex(res)
        res[ind] = x[i[ind]]
    end
    return res
end

Base.getindex(x::BioSequence, ::Colon) = copy(x)

## Setindex
Base.@propagate_inbounds function Base.setindex!(x::BioSequence, v, i::Integer)
    @boundscheck checkbounds(x, i)
    vT = convert(eltype(typeof(x)), v)
    data = encode(Alphabet(x), vT)
    encoded_setindex!(x, data, i)
end

Base.@propagate_inbounds function Base.setindex!(seq::BioSequence, x, locs::AbstractVector{<:Integer})
    @boundscheck checkbounds(seq, locs)
    @boundscheck if length(x) != length(locs)
        throw(DimensionMismatch("Attempt to assign $(length(x)) values to $(length(locs)) destinations"))
    end
    for (i, xi) in zip(locs, x)
        @boundscheck checkbounds(seq, i)
        seq[i] = xi
    end
    return seq
end

Base.@propagate_inbounds function Base.setindex!(seq::BioSequence, x, locs::AbstractVector{Bool})
    @boundscheck checkbounds(seq, locs)
    n = count(locs)
    @boundscheck if length(x) != n
        throw(DimensionMismatch("Attempt to assign $(length(x)) values to $n destinations"))
    end
    j = 0
    @inbounds for i in eachindex(locs)
        if locs[i]
            j += 1
            seq[i] = x[j]
        end
    end
    return seq
end

function Base.setindex!(seq::BioSequence, x, ::Colon)
    return setindex!(seq, x, 1:lastindex(seq))
end
