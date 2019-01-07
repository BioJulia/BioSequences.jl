# Skipmer & Kmer types
# ====================

struct Skipmer{U <: Unsigned, A <: NucleicAcidAlphabet{2}, M, N, K} <: BioSequence{A}
    bits::U
    function Skipmer{U, A, M, N, K}(x::U) where {U <: Unsigned, A <: NucleicAcidAlphabet{2}, M, N, K}
        checkskipmer(Skipmer{U, A, M, N, K})
        mask = (one(U) << (2 * K)) - 1
        return new(x & mask)
    end
end

const Kmer{U, A <: NucleicAcidAlphabet{2}, K} = Skipmer{U, A, 1, 1, K}
const DNAKmer{K} = Kmer{UInt64, DNAAlphabet{2}, K}
const RNAKmer{K} = Kmer{UInt64, RNAAlphabet{2}, K}
const BigDNAKmer{K} = Kmer{UInt128, DNAAlphabet{2}, K}
const BigRNAKmer{K} = Kmer{UInt128, RNAAlphabet{2}, K}
const DNACodon = DNAKmer{3}
const RNACodon = RNAKmer{3}

@inline encoded_data_eltype(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = U

@inline encoded_data(x::Skipmer) = x.bits

@inline cycle_len(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = N

@inline bases_per_cycle(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = M

@inline kmersize(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = K
@inline kmersize(skipmer::Skipmer) = kmersize(typeof(skipmer))

@inline capacity(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K} = div(8 * sizeof(U), 2)
@inline capacity(x::Skipmer) = capacity(typeof(x))
@inline n_unused(x::Skipmer) = capacity(x) - length(x)

_span(M, N, K) = N * (K / M - 1) + M
@inline function span(::Type{T}) where {T <: Skipmer}
    return _span(bases_per_cycle(T), cycle_len(T), kmersize(T))
end
@inline span(skipmer::T) where {T <: Skipmer} = span(typeof(skipmer))

@inline function checkskipmer(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K}
    if !(U <: Unsigned)
        throw(ArgumentError("Parameter U must be an unsigned integer type"))
    end
    if !(A <: NucleicAcidAlphabet{2})
        throw(ArgumentError("Skipmer types must have a NucleicAcidAlphabet{2} alphabet"))
    end
    capacity = div(8 * sizeof(U), 2)
    if !(1 ≤ K ≤ capacity)
        throw(ArgumentError("K must be within 1..$capacity"))
    end
    if M > N
        throw(ArgumentError("M must not be greater than N in Skipmer{A, M, N, K}"))
    end
end

include("conversion.jl")

Alphabet(::Type{Skipmer{U, A, M, N, K} where A<:NucleicAcidAlphabet{2}}) where {U, M, N, K} = Any


# Base Functions
# --------------

Base.length(x::Skipmer) = kmersize(x)

Base.summary(x::DNAKmer{k}) where {k} = string("DNA ", k, "-mer")
Base.summary(x::RNAKmer{k}) where {k} = string("RNA ", k, "-mer")
Base.summary(x::Skipmer{U, DNAAlphabet{2}, M, N, K}) where {U, M, N, K} = string("DNA Skip(", M, ", ", N, ", ", K, ")-mer")
Base.summary(x::Skipmer{U, RNAAlphabet{2}, M, N, K}) where {U, M, N, K} = string("RNA Skip(", M, ", ", N, ", ", K, ")-mer")

function Base.typemin(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K}
    checkskipmer(Skipmer{U, A, M, N, K})
    return Skipmer{U, A, M, N, K}(typemin(U))
end

function Base.typemax(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K}
    checkskipmer(Skipmer{U, A, M, N, K})
    return Skipmer{U, A, M, N, K}(typemax(U) >> (8 * sizeof(U) - 2K))
end 

function Base.rand(::Type{T}) where {T <: Skipmer}
    return T(rand(encoded_data_eltype(T)))
end

function Base.rand(::Type{T}, size::Integer) where {T <: Skipmer}
    return [rand(T) for _ in 1:size]
end

include("indexing.jl")
include("predicates.jl")
include("operations.jl")
include("transformations.jl")

Base.:-(x::Skipmer, y::Integer) = typeof(x)(encoded_data(x) - y % encoded_data_eltype(x))
Base.:+(x::Skipmer, y::Integer) = typeof(x)(encoded_data(x) + y % encoded_data_eltype(x))
Base.:+(x::Integer, y::Skipmer) = y + x

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
