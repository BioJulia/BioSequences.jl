# Skipmer & Kmer types
# ====================

primitive type Skipmer{T <: NucleicAcid, M, N, K} <: ShortSequence{64} 64 end
const Kmer{T <: NucleicAcid, K} = Skipmer{T, 1, 1, K}
const DNAKmer{K} = Kmer{DNA, K}
const RNAKmer{K} = Kmer{RNA, K}

primitive type BigSkipmer{T <: NucleicAcid, M, N, K} <: ShortSequence{128} 128 end
const BigKmer{T <: NucleicAcid, K} = BigSkipmer{T, 1, 1, K}
const BigDNAKmer{K} = BigKmer{DNA, K}
const BigRNAKmer{K} = BigKmer{RNA, K}

const DNACodon = DNAKmer{3}
const RNACodon = RNAKmer{3}

cycle_len(::Type{Skipmer{T, M, N, K}}) where {T <: NucleicAcid, M, N, K} = N
cycle_len(::Type{BigSkipmer{T, M, N, K}}) where {T <: NucleicAcid, M, N, K} = N

bases_per_cycle(::Type{Skipmer{T, M, N, K}}) where {T <: NucleicAcid, M, N, K} = M
bases_per_cycle(::Type{BigSkipmer{T, M, N, K}}) where {T <: NucleicAcid, M, N, K} = M

kmersize(::Type{Skipmer{T, M, N, K}}) where {T, M, N, K} = K
kmersize(::Type{BigSkipmer{T, M, N, K}}) where {T, M, N, K} = K
kmersize(skipmer::T) where {T <: Union{Skipmer, BigSkipmer}} = kmersize(typeof(skipmer))

_span(M, N, K) = N * (K / M - 1) + M
@inline function span(::Type{T}) where {T <: Union{Skipmer, BigSkipmer}}
    return _span(bases_per_cycle(T), cycle_len(T), kmersize(T))
end
span(skipmer::T) where T <: Union{Skipmer, BigSkipmer} = span(typeof(skipmer))

@inline function checkskipmer(::Type{Skipmer{T, M, N, K}}) where {T, M, N, K}
    if !(T <: NucleicAcid)
        throw(ArgumentError("T must be DNA or RNA in Skipmer{T, M, N, K}"))
    end
    if !(1 ≤ K ≤ 32)
        throw(ArgumentError("K must be within 1..32 in Skipmer{T, M, N, K}"))
    end
    if M > N
        throw(ArgumentError("M must not be greater than N in Skipmer{T, M, N, K}"))
    end
end

@inline function checkskipmer(::Type{BigSkipmer{T, M, N, K}}) where {T, M, N, K}
    if !(T <: NucleicAcid)
        throw(ArgumentError("T must be DNA or RNA in BigSkipmer{T, M, N, K}"))
    end
    if !(1 ≤ K ≤ 64)
        throw(ArgumentError("K must be within 1..64 in BigSkipmer{T, M, N, K}"))
    end
    if M > N
        throw(ArgumentError("M must not be greater than N in BigSkipmer{T, M, N, K}"))
    end
end 

include("conversion.jl")

BioSymbols.alphabet(::Type{T}) where {T <: Union{Skipmer{DNA}, BigSkipmer{DNA}}} = (DNA_A, DNA_C, DNA_G, DNA_T)
BioSymbols.alphabet(::Type{T}) where {T <: Union{Skipmer{RNA}, BigSkipmer{RNA}}} = (RNA_A, RNA_C, RNA_G, RNA_U)
Alphabet(::Type{Skipmer{T, M, N, K} where T<:NucleicAcid}) where {M, N, K} = Any
Alphabet(::Type{BigSkipmer{T, M, N, K} where T<:NucleicAcid}) where {M, N, K} = Any
Alphabet(::Type{T}) where {T <: Union{Skipmer{DNA}, BigSkipmer{DNA}}} = DNAAlphabet{2}()
Alphabet(::Type{T}) where {T <: Union{Skipmer{RNA}, BigSkipmer{RNA}}} = RNAAlphabet{2}()

# Base Functions
# --------------

Base.length(x::T) where {T <: Union{Skipmer, BigSkipmer}} = kmersize(x)

Base.summary(x::DNAKmer{k}) where {k} = string("DNA ", k, "-mer")
Base.summary(x::RNAKmer{k}) where {k} = string("RNA ", k, "-mer")
Base.summary(x::Skipmer{DNA, M, N, K}) where {M, N, K} = string("DNA Skip(", M, ", ", N, ", ", K, ")-mer")
Base.summary(x::Skipmer{RNA, M, N, K}) where {M, N, K} = string("RNA Skip(", M, ", ", N, ", ", K, ")-mer")

function Base.typemin(::Type{Skipmer{T, M, N, K}}) where {T, M, N, K}
    checkskipmer(Skipmer{T, M, N, K})
    return reinterpret(Skipmer{T, M, N, K}, UInt64(0))
end

function Base.typemin(::Type{BigSkipmer{T, M, N, K}}) where {T, M, N, K}
    checkskipmer(BigSkipmer{T, M, N, K})
    return reinterpret(BigSkipmer{T, M, N, K}, UInt128(0))
end

function Base.typemax(::Type{Skipmer{T, M, N, K}}) where {T, M, N, K}
    checkskipmer(Skipmer{T, M, N, K})
    return reinterpret(Skipmer{T, M, N, K}, ~UInt64(0) >> (64 - 2K))
end 

function Base.typemax(::Type{BigSkipmer{T, M, N, K}}) where {T, M, N, K}
    checkskipmer(BigSkipmer{T, M, N, K})
    return reinterpret(BigSkipmer{T, M, N, K}, ~UInt128(0) >> (128 - 2K))
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
struct KmerNeighborIterator{T, K}
    x::Kmer{T, K}
end

"""
    neighbors(kmer::Kmer)

Return an iterator through k-mers neighboring `kmer` on a de Bruijn graph.
"""
neighbors(x::Kmer{T,K}) where {T,K} = KmerNeighborIterator{T,K}(x)

Base.length(::KmerNeighborIterator) = 4
Base.eltype(::Type{KmerNeighborIterator{T,k}}) where {T,k} = Kmer{T,k}

function Base.iterate(it::KmerNeighborIterator{T, K}, i::UInt64=UInt64(0)) where {T,K}
    if i == 4
        return nothing
    else
        return Kmer{T,K}((UInt64(it.x) << 2) | i), i + 1
    end
end


# String literal
# --------------

macro kmer_str(seq)
    return DNAKmer(remove_newlines(seq))
end
