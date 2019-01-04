# Skipmer & Kmer types
# ====================

primitive type Skipmer{A <: NucleicAcidAlphabet{2}, M, N, K} <: ShortSequence{64, A} 64 end
const Kmer{A <: NucleicAcidAlphabet{2}, K} = Skipmer{A, 1, 1, K}
const DNAKmer{K} = Kmer{DNAAlphabet{2}, K}
const RNAKmer{K} = Kmer{RNAAlphabet{2}, K}

primitive type BigSkipmer{A <: NucleicAcidAlphabet{2}, M, N, K} <: ShortSequence{128, A} 128 end
const BigKmer{A <: NucleicAcidAlphabet{2}, K} = BigSkipmer{A, 1, 1, K}
const BigDNAKmer{K} = BigKmer{DNAAlphabet{2}, K}
const BigRNAKmer{K} = BigKmer{RNAAlphabet{2}, K}

const DNACodon = DNAKmer{3}
const RNACodon = RNAKmer{3}

cycle_len(::Type{Skipmer{A, M, N, K}}) where {A <: NucleicAcidAlphabet{2}, M, N, K} = N
cycle_len(::Type{BigSkipmer{A, M, N, K}}) where {A <: NucleicAcidAlphabet{2}, M, N, K} = N

bases_per_cycle(::Type{Skipmer{A, M, N, K}}) where {A <: NucleicAcidAlphabet{2}, M, N, K} = M
bases_per_cycle(::Type{BigSkipmer{A, M, N, K}}) where {A <: NucleicAcidAlphabet{2}, M, N, K} = M

kmersize(::Type{Skipmer{A, M, N, K}}) where {A, M, N, K} = K
kmersize(::Type{BigSkipmer{A, M, N, K}}) where {A, M, N, K} = K
kmersize(skipmer::T) where {T <: Union{Skipmer, BigSkipmer}} = kmersize(typeof(skipmer))

_span(M, N, K) = N * (K / M - 1) + M
@inline function span(::Type{T}) where {T <: Union{Skipmer, BigSkipmer}}
    return _span(bases_per_cycle(T), cycle_len(T), kmersize(T))
end
span(skipmer::T) where T <: Union{Skipmer, BigSkipmer} = span(typeof(skipmer))

@inline function checkskipmer(::Type{Skipmer{A, M, N, K}}) where {A, M, N, K}
    if !(A <: NucleicAcidAlphabet{2})
        throw(ArgumentError("Skipmer must have a NucleicAcidAlphabet{2}"))
    end
    if !(1 ≤ K ≤ 32)
        throw(ArgumentError("K must be within 1..32 in Skipmer{A, M, N, K}"))
    end
    if M > N
        throw(ArgumentError("M must not be greater than N in Skipmer{A, M, N, K}"))
    end
end

@inline function checkskipmer(::Type{BigSkipmer{A, M, N, K}}) where {A, M, N, K}
    if !(A <: NucleicAcidAlphabet{2})
        throw(ArgumentError("BigSkipmer must have a NucleicAcidAlphabet{2}"))
    end
    if !(1 ≤ K ≤ 64)
        throw(ArgumentError("K must be within 1..64 in BigSkipmer{T, M, N, K}"))
    end
    if M > N
        throw(ArgumentError("M must not be greater than N in BigSkipmer{T, M, N, K}"))
    end
end 

include("conversion.jl")

Alphabet(::Type{Skipmer{A, M, N, K} where A<:NucleicAcidAlphabet{2}}) where {M, N, K} = Any
Alphabet(::Type{BigSkipmer{A, M, N, K} where A<:NucleicAcidAlphabet{2}}) where {M, N, K} = Any


# Base Functions
# --------------

Base.length(x::T) where {T <: Union{Skipmer, BigSkipmer}} = kmersize(x)

Base.summary(x::DNAKmer{k}) where {k} = string("DNA ", k, "-mer")
Base.summary(x::RNAKmer{k}) where {k} = string("RNA ", k, "-mer")
Base.summary(x::Skipmer{DNAAlphabet{2}, M, N, K}) where {M, N, K} = string("DNA Skip(", M, ", ", N, ", ", K, ")-mer")
Base.summary(x::Skipmer{RNAAlphabet{2}, M, N, K}) where {M, N, K} = string("RNA Skip(", M, ", ", N, ", ", K, ")-mer")

function Base.typemin(::Type{Skipmer{A, M, N, K}}) where {A, M, N, K}
    checkskipmer(Skipmer{A, M, N, K})
    return reinterpret(Skipmer{A, M, N, K}, UInt64(0))
end

function Base.typemin(::Type{BigSkipmer{A, M, N, K}}) where {A, M, N, K}
    checkskipmer(BigSkipmer{A, M, N, K})
    return reinterpret(BigSkipmer{A, M, N, K}, UInt128(0))
end

function Base.typemax(::Type{Skipmer{A, M, N, K}}) where {A, M, N, K}
    checkskipmer(Skipmer{A, M, N, K})
    return reinterpret(Skipmer{A, M, N, K}, ~UInt64(0) >> (64 - 2K))
end 

function Base.typemax(::Type{BigSkipmer{A, M, N, K}}) where {A, M, N, K}
    checkskipmer(BigSkipmer{A, M, N, K})
    return reinterpret(BigSkipmer{A, M, N, K}, ~UInt128(0) >> (128 - 2K))
end 

function Base.rand(::Type{T}) where {T <: Union{Skipmer, BigSkipmer}}
    return T(rand(encoded_data_eltype(T)))
end

function Base.rand(::Type{T}, size::Integer) where {T <: Union{Skipmer, BigSkipmer}}
    return [rand(T) for _ in 1:size]
end


# K-mer neighbor
# --------------

# neighbors on a de Bruijn graph
struct KmerNeighborIterator{A, K}
    x::Kmer{A, K}
end

"""
    neighbors(kmer::Kmer)

Return an iterator through k-mers neighboring `kmer` on a de Bruijn graph.
"""
neighbors(x::Kmer{A,K}) where {A,K} = KmerNeighborIterator{A,K}(x)

Base.length(::KmerNeighborIterator) = 4
Base.eltype(::Type{KmerNeighborIterator{A,K}}) where {A,K} = Kmer{A,K}

function Base.iterate(it::KmerNeighborIterator{A, K}, i::UInt64=UInt64(0)) where {A,K}
    if i == 4
        return nothing
    else
        return Kmer{A,K}((UInt64(it.x) << 2) | i), i + 1
    end
end


# String literal
# --------------

macro kmer_str(seq)
    return DNAKmer(remove_newlines(seq))
end
