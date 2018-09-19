# BioSequence
# ===========
#
# Abstract biological sequence type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


# Type Definition
# ---------------

"""
# An abstract biological sequence type.

## Required interface

Any subtype of `BioSequence` should implement the following methods:

* `encoded_data(seq)`: return the thing containing the encoded elements of `seq`. 
* `Alphabet(typeof(seq))`: return the type of the alphabet of the type of `seq`.
* `Base.length(seq)`: return the length of `seq`.
* `inbounds_getindex(seq, i)`: return the i'th element of `seq`, without checking bounds.
* `find_next_ambiguous(seq)`: return the index of the next ambiguous symbol in the
  sequence.

  Any sequence type that adheres to the above can be expected to have the
  following implemented "for free".

## Provided, BioSequences traits

* `encoded_data_eltype(seq)`: return the element type of the thing contaning encoded elements of `seq`.
* `Alphabet(seq)`: return the type of the alphabet of `seq`.
* `BitsPerSymbol(seq)`: return the number of bits `seq` uses to encode each element, as a type trait.
* `bits_per_symbol(seq)`: return the number of bits `seq` uses to encode each element, as a value.
* `symbols_per_data_element(seq)`: return the number of elements that can be encoded into each element in `encoded_data(seq)`.
* `bitindex(seq, i)`: return the BitIndex value corresponding to the i'th element in `seq`.
* `firstbitindex(seq)`: return the BitIndex value corresponding to the first element in `seq`.
* `lastbitindex(seq)`: return the BitIndex value corresponding to the last element in `seq`.

## Provided Base traits & methods

* `Base.eltype(typeof(seq))`: return the element type for the type of `seq`.
* `Base.eltype(seq)`: return the element type of `seq`.
* `Base.size(seq)`: Return the size of `seq`.
* `Base.firstindex(seq)`: Return the first index of `seq`.
* `Base.lastindex(seq::S)`: Return the last index of `seq`.
* `Base.eachindex(seq::S)`: Return an iterator over all indicies of `seq`.
* `Base.getindex(seq::S, i::Integer)`: Return the biological symbol `seq`
  contains at position `i`.
* `Base.:(==)(seq1::S, seq2::S)`: Check if `seq1` and `seq2` can be considered
  equal in value.
* `Base.isless(seq1::S, seq2::S)`: Check if `seq1` is considered lesser in
  value than `seq2`.
* `Base.cmp(seq1::S, seq2::S)`: Returns -1, 0, or 1 depending on whether
  `seq1` is less than, equal to, or greater than `seq2`, respectively.
* `Base.findnext(seq::S, val, start::Integer)`: Return the index of the next
  occurance of `val` in `seq`, starting from index `start`.
* `Base.findprev(seq::S, val, start::Integer)`: Return the index of the
  previous occurance of `val` in `seq`, starting from index `start`.
* `Base.findfirst(seq::S, val)`: Return the index of the first occurance of
  `val` in `seq`.
* `Base.findlast(seq::S, val)`: Return the index of the final occurance of
  `val` in `seq`.
* `Base.isempty(seq::S)`: Determine whether `seq` has no elements.
"""
abstract type BioSequence end

# This is useful for obscure reasons. We use SeqRecord{BioSequence} for reading
# sequence in an undetermined alphabet, but a consequence that we need to be
# able to construct a `Sequence`.
function BioSequence()
    return DNASequence()
end

# Required traits and methods
# ---------------------------

"""
Return the data member of `seq` that stores the encoded sequence data.
"""
@inline function encoded_data(seq::BioSequence)
    error(
        string(
            "encoded_data has not been defined for BioSequence type: ",
            typeof(seq),
            ". It is required for any BioSequence subtype."
        )
    )
end

encoded_data_eltype(seq::BioSequence) = eltype(encoded_data(seq))

"""
Return the `Alpahbet` type defining the possible biological symbols
and their encoding for a given biological sequence.
"""
@inline function Alphabet(::Type{S}) where S <: BioSequence
    error(string("This sequence type trait has not been defined for BioSequence type: ", S))
end

# Provided traits and methods
# ---------------------------

# This version of Alphabet is automatically defined for any BioSequence type.
# Is more for conveinience.
@inline function Alphabet(seq::BioSequence)
    return Alphabet(typeof(seq))
end

BitsPerSymbol(seq::BioSequence) = BitsPerSymbol(Alphabet(seq))
bits_per_symbol(seq::BioSequence) = bits_per_symbol(Alphabet(seq))

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


# Provided Base methods
# ---------------------

Base.eltype(::Type{T}) where T <: BioSequence = eltype(Alphabet(T))
Base.eltype(seq::BioSequence) = eltype(Alphabet(seq))
Base.size(seq::BioSequence) = (length(seq),)

include("indexing.jl")
include("conversion.jl")
include("predicates.jl")

# Finders
# -------

function Base.findnext(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:lastindex(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findprev(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

Base.findfirst(val, seq::BioSequence) = findnext(val, seq, 1)
Base.findlast(val, seq::BioSequence) = findprev(val, seq, lastindex(seq))

include("printing.jl")
include("operations.jl")
include("transformations.jl")
