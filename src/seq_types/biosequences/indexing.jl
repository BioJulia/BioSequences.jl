# Indexing
# ========
#
# Indexing methods for mutable biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

Base.eltype(::Type{T}) where T <: BioSequence = eltype(alphabet_t(T))
Base.eltype(seq::BioSequence) = eltype(alphabet_t(seq))
Base.size(seq::BioSequence) = (length(seq),)
Base.endof(seq::BioSequence) = length(seq)
Base.eachindex(seq::BioSequence) = 1:endof(seq)


# Bounds checking
# ---------------

@inline function Base.checkbounds(seq::BioSequence, i::Integer)
    if 1 ≤ i ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, i))
end

@inline function Base.checkbounds(seq::BioSequence, locs::AbstractVector{Bool})
    if length(seq) == length(locs)
        return true
    end
    throw(BoundsError(seq, locs))
end

@inline function Base.checkbounds(seq::BioSequence, locs::AbstractVector)
    for i in locs
        checkbounds(seq, i)
    end
    return true
end

@inline function Base.checkbounds(seq::BioSequence, range::UnitRange)
    if 1 ≤ range.start && range.stop ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, range))
end

function checkdimension(from::Integer, to::Integer)
    if from == to
        return true
    end
    throw(
        DimensionMismatch(
            string("attempt to assign ", from, " elements to ", to, " elements")
        )
    )
end

function checkdimension(seq::BioSequence, locs::AbstractVector)
    return checkdimension(length(seq), length(locs))
end

function checkdimension(seq::BioSequence, locs::AbstractVector{Bool})
    return checkdimension(length(seq), sum(locs))
end

@inline function inbounds_getindex(seq::BioSequence, i::Integer)
    bidx = bitindex(seq, i)
    encoded_symbol = extract_encoded_symbol(bidx, encoded_data(seq))
    return decode(alphabet_t(seq), encoded_symbol)
end

@inline function Base.getindex(seq::BioSequence, i::Integer)
    @boundscheck checkbounds(seq, i)
    return inbounds_getindex(seq, i)
end
