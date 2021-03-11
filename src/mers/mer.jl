###
### A short immutable sequence type for representing Kmers and Skipmers.
###

abstract type AbstractMer{A<:NucleicAcidAlphabet{2},K} <: BioSequence{A} end

primitive type Mer{A<:NucleicAcidAlphabet{2},K} <: AbstractMer{A,K} 64 end
primitive type BigMer{A<:NucleicAcidAlphabet{2},K} <: AbstractMer{A,K} 128 end

# Aliases
const DNAMer{K} = Mer{DNAAlphabet{2},K}
const RNAMer{K} = Mer{RNAAlphabet{2},K}
const DNAKmer   = DNAMer{31}
const RNAKmer   = RNAMer{31}

const BigDNAMer{K} = BigMer{DNAAlphabet{2},K}
const BigRNAMer{K} = BigMer{RNAAlphabet{2},K}
const BigDNAKmer   = BigMer{DNAAlphabet{2},63}
const BigRNAKmer   = BigMer{RNAAlphabet{2},63}

const DNACodon = DNAMer{3}
const RNACodon = RNAMer{3}

###
### Base Functions
###

@inline encoded_data_type(::Type{Mer{A,K}}) where {A,K} = UInt64
@inline encoded_data_type(::Type{BigMer{A,K}}) where {A,K} = UInt128
@inline encoded_data_type(x::AbstractMer) = encoded_data_type(typeof(x))
@inline encoded_data(x::AbstractMer) = reinterpret(encoded_data_type(typeof(x)), x)
@inline capacity(::Type{U}) where {U<:Unsigned} = div(8 * sizeof(U), 2)
@inline capacity(mer::Type{T}) where {T<:AbstractMer} = capacity(encoded_data_type(mer))
@inline capacity(mer::AbstractMer) = capacity(typeof(mer))
@inline ksize(::Type{T}) where {A,K,T<:AbstractMer{A,K}} = K
@inline n_unused(::Type{T}) where {T<:AbstractMer} = capacity(T) - ksize(T)
@inline n_unused(x::AbstractMer) = n_unused(typeof(x))

@inline function checkmer(::Type{T}) where {T<:AbstractMer}
    c = capacity(T)
    if !(1 ≤ ksize(T) ≤ c)
        throw(ArgumentError("K must be within 1..$c"))
    end
end

"""
    masked(::Type{<:AbstractMer}, x::Unsigned)

Re-interpret `x` as a `T`, setting any uncoding bits to zero.
Throws an error if `T` is un-instantiable, e.g. `DNAMer{-1}`.
"""
function masked(::Type{T}, x::UInt64) where {T <: Mer}
    checkmer(T)
    reinterpret(T, x & ((one(UInt64) << (2 * ksize(T))) - 1))
end

function masked(::Type{T}, x::UInt128) where {T <: BigMer}
    checkmer(T)
    reinterpret(T, x & ((one(UInt128) << (2 * ksize(T))) - 1))
end

@inline Base.length(x::AbstractMer{A,K}) where {A,K} = ksize(typeof(x))
@inline Base.unsigned(x::AbstractMer) = encoded_data(x)
@inline Base.summary(::AbstractMer{DNAAlphabet{2},K}) where {K} = string("DNA ", K, "-mer")
@inline Base.summary(::AbstractMer{RNAAlphabet{2},K}) where {K} = string("RNA ", K, "-mer")

function Base.typemin(::Type{T}) where {T<:AbstractMer}
    checkmer(T)
    return reinterpret(T, typemin(encoded_data_type(T)))
end

function Base.typemax(::Type{T}) where {T<:AbstractMer}
    checkmer(T)
    U = encoded_data_type(T)
    return reinterpret(T, typemax(U) >> (8 * sizeof(U) - 2ksize(T)))
end 

function Base.rand(::Type{T}) where {T<:AbstractMer}
    return masked(T, rand(encoded_data_type(T)))
end

function Base.rand(::Type{T}, size::Integer) where {T<:AbstractMer}
    return [rand(T) for _ in 1:size]
end

Alphabet(::Type{Mer{A,K} where A<:NucleicAcidAlphabet{2}}) where {K} = Any

function Mer{A,K}(x::UInt64) where {A,K}
    checkmer(Mer{A,K})
    maxval = (one(UInt64) << (2 * K)) - 1
    x > maxval && throw(DomainError(x, "Mer with K=$K is too small for $x"))
    return reinterpret(Mer{A,K}, x)
end

function BigMer{A,K}(x::UInt128) where {A,K}
    checkmer(BigMer{A,K})
    maxval = (one(UInt128) << (2 * K)) - 1
    x > maxval && throw(DomainError(x, "Mer with K=$K is too small for $x"))
    return reinterpret(BigMer{A,K}, x)
end

###
### Conversion
###

include("conversion.jl")
include("indexing.jl")
include("predicates.jl")
include("counting.jl")
include("transformations.jl")

###
### Mer de-bruijn neighbors
###

# Neighbors on a de Bruijn graph
struct MerNeighborIterator{S<:AbstractMer}
    x::S
end

"""
    neighbors(skipmer::S) where {S<:AbstractMer}

Return an iterator through skip-mers neighboring `skipmer` on a de Bruijn graph.
"""
neighbors(skipmer::AbstractMer) = MerNeighborIterator{typeof(skipmer)}(skipmer)

Base.length(::MerNeighborIterator) = 4
Base.eltype(::Type{MerNeighborIterator{S}}) where {S<:AbstractMer} = S

function Base.iterate(it::MerNeighborIterator{S}, i = zero(encoded_data_type(S))) where {S<:AbstractMer}
    if i == 4
        return nothing
    else
        return masked(S, (encoded_data(it.x) << 2) | i), i + 1
    end
end

###
### String literals
###

macro mer_str(seq, flag)
    seq′ = remove_newlines(seq)
    if flag == "dna" || flag == "d"
        return DNAMer(seq′)
    elseif flag == "rna" || flag == "r"
        return RNAMer(seq′)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

macro mer_str(seq)
    return DNAMer(remove_newlines(seq))
end

macro bigmer_str(seq, flag)
    seq′ = remove_newlines(seq)
    if flag == "dna" || flag == "d"
        return BigDNAMer(seq′)
    elseif flag == "rna" || flag == "r"
        return BigRNAMer(seq′)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

macro bigmer_str(seq)
    return BigDNAMer(remove_newlines(seq))
end
