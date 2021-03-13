###
### An abstract biological sequence type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    BioSequence{A <: Alphabet}

`BioSequence` is the main abstract type of `BioSequences`.
Its subtypes are characterized by:
* Being a linear container type with random access and indices `Base.OneTo(length(x))`.
* Containing zero or more coding elements of type `encoded_data_eltype(typeof(x))`.
* Being associated with an `Alphabet`, `A` by being a subtype of `BioSequence{A}`.

A `BioSequence{A}` is indexed by an integer. The biosequence subtype, the index
and the alphabet `A` determine how to extract the internal encoding data.
The alphabet decides how to decode the data to the element type of the biosequence.
Hence, the element type and container type of a `BioSequence` are separated.

Subtypes `T` of `BioSequence` must implement the following, with `E` begin an
encoded data type:

* `Base.length(::T)::Int`
* `encoded_data_eltype(::Type{T})::Type{E}`
* `extract_encoded_element(::T, ::Integer)::E`
* T must be able to be constructed from any iterable with `length` defined and
with a known, compatible element type.

Furthermore, mutable sequences should implement
* `encoded_setindex!(::T, ::E, ::Integer)`
* `T(undef, ::Int)`
* resize!(::T, ::Int)
* empty(::Type{T})

For compatibility with existing `Alphabet`s, the encoded data eltype must be `UInt`.
"""
abstract type BioSequence{A<:Alphabet} end

Base.eachindex(x::BioSequence) = Base.OneTo(length(x))
Base.firstindex(::BioSequence) = 1
Base.lastindex(x::BioSequence) = length(x)
Base.keys(seq::BioSequence) = eachindex(seq)
Base.nextind(::BioSequence, i::Integer) = Int(i) + 1
Base.prevind(::BioSequence, i::Integer) = Int(i) - 1
Base.size(x::BioSequence) = (length(x),)
Base.eltype(::Type{<:BioSequence{A}}) where {A <: Alphabet} = eltype(A)
Alphabet(::Type{<:BioSequence{A}}) where {A <: Alphabet} = A()
Alphabet(x::BioSequence) = Alphabet(typeof(x))

# Specific biosequences
const NucleotideSeq = BioSequence{<:NucleicAcidAlphabet}
const AminoAcidSeq = BioSequence{AminoAcidAlphabet}

# The generic functions for any BioSequence...
include("indexing.jl")
include("conversion.jl")
include("predicates.jl")
include("find.jl")
include("printing.jl")
include("transformations.jl")
include("counting.jl")
