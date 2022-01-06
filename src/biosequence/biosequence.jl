###
### An abstract biological sequence type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    BioSequence{A <: Alphabet}

`BioSequence` is the main abstract type of `BioSequences`.
It abstracts over the internal representation of different biological sequences,
and is parameterized by an `Alphabet`, which controls the element type.

# Extended help
Its subtypes are characterized by:
* Being a linear container type with random access and indices `Base.OneTo(length(x))`.
* Containing zero or more internal data elements of type `encoded_data_eltype(typeof(x))`.
* Being associated with an `Alphabet`, `A` by being a subtype of `BioSequence{A}`.

A `BioSequence{A}` is indexed by an integer. The biosequence subtype, the index
and the alphabet `A` determine how to extract the internal encoded data.
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
* `resize!(::T, ::Int)`

For compatibility with existing `Alphabet`s, the encoded data eltype must be `UInt`.
"""
abstract type BioSequence{A<:Alphabet} end

"""
    has_interface(::Type{BioSequence}, ::T, syms::Vector, mutable::Bool, compat::Bool=true)

Check if type `T` conforms to the `BioSequence` interface. A `T` is constructed from the vector
of element types `syms` which must not be empty.
If the `mutable` flag is set, also check the mutable interface.
If the `compat` flag is set, check for compatibility with existing alphabets.
"""
function has_interface(
    ::Type{BioSequence},
    ::Type{T},
    syms::Vector,
    mutable::Bool,
    compat::Bool=true
) where {T <: BioSequence}
    try
        isempty(syms) && error("Vector syms must not be empty")
        first(syms) isa eltype(T) || error("Vector is of wrong element type")
        seq = T(syms)
        length(seq) > 0 || return false
        E = encoded_data_eltype(T)
        e = extract_encoded_element(seq, 1)
        e isa E || return false
        (!compat || E == UInt) || return false
        copy(seq) isa typeof(seq) || return false
        if mutable
            encoded_setindex!(seq, e, 1)
            T(undef, 5) isa T || return false
            isempty(resize!(seq, 0)) || return false
        end
    catch error
        error isa MethodError && return false
        rethrow(error)
    end
    return true
end

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
BitsPerSymbol(x::BioSequence) = BitsPerSymbol(Alphabet(typeof(x)))

function Base.similar(seq::BioSequence, len::Integer=length(seq))
    return typeof(seq)(undef, len)
end

# Fast path for iterables we know are stateless
function join!(seq::BioSequence, it::Union{Vector, Tuple, Set})
    _join!(resize!(seq, sum(length, it, init=0)), it, Val(true))
end

"""
    join!(seq::BioSequence, iter)

Concatenate all biosequences in `iter` into `seq`, resizing it to fit.

# Examples
```
julia> join(LongDNA(), [dna"TAG", dna"AAC"])
6nt DNA Sequence:
TAGAAC
```

see also [`join`](@ref)
"""
join!(seq::BioSequence, it) = _join!(seq, it, Val(false))

# B is whether the size of the destination seq is already
# known to be the final size
function _join!(seq::BioSequence, it, ::Val{B}) where B
    index = 1
    for i in it
        B || resize!(seq, length(seq) + length(i))
        copyto!(seq, index, i, 1, length(i))
        index += length(i)
    end
    seq
end

"""
    join(::Type{T <: BioSequence}, seqs)

Concatenate all the `seqs` to a biosequence of type `T`.

# Examples
```
julia> join(LongDNA, [dna"TAG", dna"AAC"])
6nt DNA Sequence:
TAGAAC
```

see also [`join!`](@ref)
"""
function Base.join(::Type{T}, it::Union{Vector, Tuple, Set}) where {T <: BioSequence}
    _join!(T(undef, sum(length, it, init=0)), it, Val(true))
end

function Base.join(::Type{T}, it) where {T <: BioSequence}
    _join!(empty(T), it, Val(false))
end

Base.repeat(chunk::BioSequence, n::Integer) = join(typeof(chunk), (chunk for i in 1:n))
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
"""
An alias for `BioSequence{<:NucleicAcidAlphabet}`
"""
const NucleotideSeq = BioSequence{<:NucleicAcidAlphabet}

"An alias for `BioSequence{<:NucleicAcidAlphabet}`"
const NucSeq{N} = BioSequence{<:NucleicAcidAlphabet{N}}

"An alias for `BioSequence{DNAAlphabet{N}}`"
const DNASeq{N} = BioSequence{DNAAlphabet{N}}

"An alias for `BioSequence{RNAAlphabet{N}}`"
const RNASeq{N} = BioSequence{RNAAlphabet{N}}

"""
An alias for `BioSequence{AminoAcidAlphabet}`
"""
const AASeq = BioSequence{AminoAcidAlphabet}

# The generic functions for any BioSequence...
include("indexing.jl")
include("conversion.jl")
include("predicates.jl")
include("find.jl")
include("printing.jl")
include("transformations.jl")
include("counting.jl")
include("copying.jl")
