###
### Indexing
###
###
### Indexing methods for mutable biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

#=
Indexing into BioSequences works as follows:

A `Base.getindex` method is provided for `BioSequence`.
It's job is to, do boundschecking. And to then call `inbounds_getindex`.

It is the task of `inbounds_getindex` to extract from a sequence, the
appropriate bits, and pass them to an appropriate `decode` method.
=#

@inline function symbols_per_data_element(seq::BioSequence)
    return div(8 * sizeof(encoded_data_eltype(seq)), bits_per_symbol(seq))
end

@inline function bitindex(seq::BioSequence, i::Integer)
    return bitindex(BitsPerSymbol(seq), encoded_data_eltype(seq), i)
end

@inline function firstbitindex(seq::BioSequence)
    return bitindex(seq, 1)
end

@inline function lastbitindex(seq::BioSequence)
    return bitindex(seq, lastindex(seq))
end

@inline function bindata_mask(seq::BioSequence)
    return bitmask(Alphabet(seq))
end

function eachbitindex(seq::BioSequence, from = firstindex(seq), to = lastindex(seq))
    return bitindex(seq, from):bits_per_symbol(seq):bitindex(seq, to)
end


"""
    inbounds_getindex(seq::BioSequence, i::BitIndex)

An inbounds_getindex method defined by default for all BioSequences.

It is the job of this method to take a sequence and a bitindex, extract the
encoded bits and decode them into a symbol.

!!! note
    This method might be overloaded for specialised subtypes of BioSequence.
"""
@inline function inbounds_getindex(seq::BioSequence, bidx::BitIndex)
    encoded_symbol = extract_encoded_element(bidx, encoded_data(seq))
    return decode(Alphabet(seq), encoded_symbol)
end

"""
    inbounds_getindex(seq::BioSequence, i::Integer)

An inbounds_getindex method defined by default for all BioSequences.

It is the job of this method to convert the integer index to a `BitIndex`, and
call the method of inbounds_getindex, for the sequence and `BitIndex`
combination.

!!! note
    This method might be overloaded for specialised subtypes of BioSequence.
"""
@inline function inbounds_getindex(seq::BioSequence, i::Integer)
    bidx = bitindex(seq, i)
    return inbounds_getindex(seq::BioSequence, bidx)
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

@inline function unsafe_setindex!(s::BioSequence, v, i::Integer)
    return unsafe_setindex!(s, v, bitindex(s, i))
end

@inline function unsafe_setindex!(s::BioSequence, v, i::BitIndex)
    v_ = convert(eltype(s), v)
    encoded = encode(Alphabet(s), v_)
    return encoded_setindex!(s, encoded, i)
end

###
### Base methods
###

# Provided Base methods for indexing any BioSequence type. May be overloaded if
# needed but it is unlikely.

Base.eltype(::Type{T}) where T <: BioSequence = eltype(Alphabet(T))
Base.eltype(seq::BioSequence) = eltype(Alphabet(seq))
Base.size(seq::BioSequence) = (length(seq),)

Base.firstindex(seq::BioSequence) = 1
Base.lastindex(seq::BioSequence) = length(seq)
Base.eachindex(seq::BioSequence) = Base.OneTo(lastindex(seq))
Base.keys(seq::BioSequence) = eachindex(seq)
Base.nextind(::BioSequence, i::Integer) = Int(i) + 1
Base.prevind(::BioSequence, i::Integer) = Int(i) - 1

# Bounds checking...
@inline function Base.checkbounds(seq::BioSequence, i::Integer)
    if 1 ≤ i ≤ lastindex(seq)
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
    if 1 ≤ range.start && range.stop ≤ lastindex(seq)
        return true
    end
    throw(BoundsError(seq, range))
end

@inline function Base.getindex(seq::BioSequence, i::Integer)
    @boundscheck checkbounds(seq, i)
    return inbounds_getindex(seq, i)
end

@inline function Base.setindex!(s::BioSequence, v, i::Integer)
    @boundscheck checkbounds(s, i)
    return unsafe_setindex!(s, v, i)
end

function Base.setindex!(seq::BioSequence, x, locs::AbstractVector{<:Integer})
    @boundscheck checkbounds(seq, locs)
    return unsafe_setindex!(seq, x, locs)
end

function Base.setindex!(seq::BioSequence, x::BioSequence, locs::AbstractVector{<:Integer})
    @boundscheck checkbounds(seq, locs)
    checkdimension(x, locs)
    @inbounds for i in eachindex(locs)
        unsafe_setindex!(seq, x[i], locs[i])
    end
    seq
end

function Base.setindex!(seq::BioSequence, x::BioSequence, locs::AbstractVector{Bool})
    @boundscheck checkbounds(seq, locs)
    checkdimension(x, locs)
    j = 0
    @inbounds for i in eachindex(locs)
        if locs[i]
            j += 1
            unsafe_setindex!(seq, x[j], i)
        end
    end
    seq
end

function Base.setindex!(seq::BioSequence, other::BioSequence, ::Colon)
    return setindex!(seq, other, 1:lastindex(seq))
end

function Base.setindex!(seq::BioSequence, x, ::Colon)
    return setindex!(seq, x, 1:lastindex(seq))
end

@inline function Base.iterate(seq::BioSequence, i::Int = firstindex(seq))
    if i > lastindex(seq)
        return nothing
    else
        return inbounds_getindex(seq, i), i + 1
    end
end
