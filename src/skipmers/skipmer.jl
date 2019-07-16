###
### A short immutable sequence type for representing Kmers and Skipmers.
###

abstract type AbstractMer{A<:NucleicAcidAlphabet{2},K} <: BioSequence{A} end

primitive type Mer{A<:NucleicAcidAlphabet{2},K} <: AbstractMer{A,K} 64 end
primitive type BigMer{A<:NucleicAcidAlphabet{2},K} <: AbstractMer{A,K} 128 end

# Aliases
const DNAMer{K} = Mer{DNAAlphabet{2},K}
const RNAMer{K} = Mer{RNAAlphabet{2},K}
const BigDNAMer{K} = BigMer{DNAAlphabet{2},K}
const BigRNAMer{K} = BigMer{RNAAlphabet{2},K}

const DNACodon = DNAMer{3}
const RNACodon = RNAMer{3}

const DEFAULT_SMALL_K = 27
const DEFAULT_LARGE_K = 31

###
### Base Functions
###

@inline encoded_data_type(::Type{Mer{A,K}}) where {A,K} = UInt64
@inline encoded_data_type(::Type{BigMer{A,K}}) where {A,K} = UInt128
@inline encoded_data_type(x::AbstractMer) = encoded_data_type(typeof(x))
@inline encoded_data(x::AbstractMer) = reintepret(encoded_data_type(typeof(x)), x)
@inline capacity(x::Type{U}) where {U<:Unsigned} = div(8 * sizeof(U), 2)
@inline capacity(mer::Type{T}) where {T<:AbstractMer} = encoded_data_type(mer)
@inline capacity(mer::AbstractMer) = capacity(typeof(mer))

@inline K(::Type{T}) where {A,K,T<:AbstractMer{A,K}} = K
@inline Base.length(x::AbstractMer{A,K}) where {A,K} = K(typeof(x))
@inline Base.unsigned(x::AbstractMer) = encoded_data(x)
@inline Base.summary(x::AbstractMer{DNAAlphabet{2},K}) where {K} = string("DNA ", K, "-mer")
@inline Base.summary(x::AbstractMer{RNAAlphabet{2},K}) where {K} = string("RNA ", K, "-mer")
Base.typemin(::Type{T}) where {T<:AbstractMer} = T(typemin(encoded_data_type(T)))

function Base.typemax(::Type{T}) where {T<:AbstractMer}
    U = encoded_data_type(T)
    return T(typemax(U) >> (8 * sizeof(U) - 2K))
end 

function Base.rand(::Type{T}) where {T<:AbstractMer}
    return T(rand(encoded_data_type(T)))
end

function Base.rand(::Type{T}, size::Integer) where {T<:AbstractMer}
    return [rand(T) for _ in 1:size]
end

Base.:-(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) - y % encoded_data_type(x))
Base.:+(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) + y % encoded_data_type(x))
Base.:+(x::AbstractMer, y::AbstractMer) = y + x

Alphabet(::Type{Mer{A,K} where A<:NucleicAcidAlphabet{2}}) where {K} = Any

###
### Conversion
###

function Mer{A,K}(x::UInt64) where {A,K}
    mask = (one(UInt64) << (2 * K)) - 1
    return reinterpret(Mer{A,K}, x & mask)
end

function BigMer{A,K}(x::UInt128) where {A,K}
    mask = (one(UInt128) << (2 * K)) - 1
    return reinterpret(BigMer{A,K}, x & mask)
end

include("conversion.jl")
include("indexing.jl")
include("predicates.jl")
include("operations.jl")
include("transformations.jl")


#@inline encoded_data(x::Skipmer) = x.bits

#@inline cycle_len(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = N

#@inline bases_per_cycle(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = M

#@inline kmersize(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = K
#@inline kmersize(skipmer::Skipmer) = kmersize(typeof(skipmer))

#@inline skipmer_capacity(::Type{U}) where {U<:Unsigned} = div(8 * sizeof(U), 2)
#@inline capacity(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = div(8 * sizeof(U), 2)
#@inline capacity(x::Skipmer) = capacity(typeof(x))
#@inline n_unused(x::Skipmer) = capacity(x) - length(x)

#_span(M::Integer, N::Integer, K::Integer) = UInt(ceil(N * (K / M - 1) + M))
#@inline function span(::Type{T}) where {T <: Skipmer}
#    return _span(bases_per_cycle(T), cycle_len(T), kmersize(T))
#end
#@inline span(skipmer::T) where {T <: Skipmer} = span(typeof(skipmer))
#
#@inline function checkskipmer(::Type{SK}) where {SK<:Skipmer}
#    c = capacity(SK)
#    if !(1 ≤ kmersize(SK) ≤ c)
#        throw(ArgumentError("K must be within 1..$c"))
#    end
#    if bases_per_cycle(SK) > cycle_len(SK)
#        throw(ArgumentError("M must not be greater than N in Skipmer{A, M, N, K}"))
#    end
#end

#=
# K-mer neighbor
# --------------

# Neighbors on a de Bruijn graph
struct SkipmerNeighborIterator{S <: Skipmer}
    x::S
end

"""
    neighbors(skipmer::S) where {S <: Skipmer}

Return an iterator through skip-mers neighboring `skipmer` on a de Bruijn graph.
"""
neighbors(skipmer::S) where {S <: Skipmer} = SkipmerNeighborIterator{S}(skipmer)

Base.length(::SkipmerNeighborIterator) = 4
Base.eltype(::Type{SkipmerNeighborIterator{S}}) where {S <: Skipmer} = S

function Base.iterate(it::SkipmerNeighborIterator{S}, i::UInt64 = UInt64(0)) where {S <: Skipmer}
    if i == 4
        return nothing
    else
        return S((encoded_data(it.x) << 2) | i), i + 1
    end
end


# String literal
# --------------

macro kmer_str(seq)
    return DNAKmer(remove_newlines(seq))
end
=#