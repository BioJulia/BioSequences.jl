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
* `copy(::T)`
* T must be able to be constructed from any iterable with `length` defined and
with a known, compatible element type.

Furthermore, mutable sequences should implement
* `encoded_setindex!(::T, ::E, ::Integer)`
* `T(undef, ::Int)`
* resize!(::T, ::Int)

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
Base.eltype(x::BioSequence) = eltype(typeof(x))
Alphabet(::Type{<:BioSequence{A}}) where {A <: Alphabet} = A()
Alphabet(x::BioSequence) = Alphabet(typeof(x))
Base.isempty(x::BioSequence) = iszero(length(x))
Base.empty(::Type{T}) where {T <: BioSequence} = T(eltype(T)[])
Base.empty(x::BioSequence) = empty(typeof(x))

function Base.similar(seq::BioSequence, len::Integer=length(seq))
    return typeof(seq)(undef, len)
end

# Fast path for iterables we know are stateless
function join!(seq::BioSequence, it::Union{Vector, Tuple, Set})
    _join(resize!(seq, sum(length, it)), it, Val(true))
end

join!(seq::BioSequence, it) = _join!(seq, it, Val(false))

function _join!(seq::BioSequence, it, ::Val{B}) where B
    index = 1
    for i in it
        B || resize!(seq, length(seq) + length(i))
        copyto!(seq, index, i, 1, length(i))
        index += length(i)
    end
    seq
end

function Base.join(::Type{T}, it::Union{Vector, Tuple, Set}) where {T <: BioSequence}
    _join!(T(undef, sum(length, it)), it, Val(true))
end

function Base.join(::Type{T}, it) where {T <: BioSequence}
    _join!(empty(T), it, Val(false))
end

Base.repeat(chunk::BioSequence, n::Integer) = join((chunk for i in 1:n))
Base.:^(x::BioSequence, n::Integer) = repeat(x, n)

# Concatenation and Base.repeat operators
function Base.:*(chunks::BioSequence...)
    T = typeof(first(chunks))
    join(T, chunks)
end

"""
    encoded_data_eltype(::Type{<:BioSequence})

Returns the element type of the encoded data of the `BioSequence`.
This is the return type of `extract_encoded_element`, i.e. the data
type that stores the biological symbols in the biosequence.

See also: [`BioSequence`](@ref) 
"""
function encoded_data_eltype end

"""
    extract_encoded_element(::BioSequence{A}, i::Integer)

Returns the encoded element at position `i`. This data can be
decoded using `decode(A(), data)` to yield the element type of
the biosequence.

See also: [`BioSequence`](@ref) 
"""
function extract_encoded_element end


"""
    encoded_setindex!(seq::BioSequence, x::E, i::Integer)

Given encoded data `x` of type `encoded_data_eltype(typeof(seq))`,
sets the internal sequence data at the given index.

See also: [`BioSequence`](@ref) 
"""
function encoded_setindex! end

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
include("copying.jl")
