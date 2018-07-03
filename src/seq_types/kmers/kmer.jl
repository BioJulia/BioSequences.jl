# K-mer
# =====
#
# Compact k-mer sequence type.
#
# A Kmer is a sequence ≤ 32nt, without any 'N's, packed in a single 64 bit
# value.  While `GeneralSequence` is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Representation
# --------------
#
# Four kinds of nucleotides are encoded as follows:
#
#   nucleotide | binary
#   ---------- | ------
#       A      |   00
#       C      |   01
#       G      |   10
#     T / U    |   11
#
# NucleicAcids are filled from MSBs to LSBs and right-aligned so that all k-mers
# are lexicographically ordered. For example, the memory layout of "TACG" is:
#   64-bit: 0b 00 00 … 00 11 00 01 10
#    4-mer:                T  A  C  G

primitive type Kmer{T <: NucleicAcid, K} <: BioSequence 64 end

const DNAKmer{K} = Kmer{DNA, K}
const RNAKmer{K} = Kmer{RNA, K}
const DNACodon = DNAKmer{3}
const RNACodon = RNAKmer{3}


function Kmer(nts::T...) where {T<:NucleicAcid}
    return make_kmer(nts)
end

function DNACodon(x::DNA, y::DNA, z::DNA)
    return make_kmer((x, y, z))
end

function RNACodon(x::RNA, y::RNA, z::RNA)
    return make_kmer((x, y, z))
end


# Conversion
# ----------

function Kmer{T, K}(x::UInt64) where {T, K}
    checkkmer(Kmer{T,K})
    mask = ~UInt64(0) >> (64 - 2K)
    return reinterpret(Kmer{T, K}, x & mask)
end

UInt64(x::Kmer) = reinterpret(UInt64, x)
Base.convert(::Type{UInt64}, x::Kmer) = reinterpret(UInt64, x)

DNAKmer{K}(x::RNAKmer{K}) where {K} = reinterpret(DNAKmer{K}, x)
RNAKmer{K}(x::DNAKmer{K}) where {K} = reinterpret(RNAKmer{K}, x)

function Kmer{T,K}(seq::AbstractString) where {T,K}
    return make_kmer(Kmer{T,K}, seq)
end

function Kmer{T,K}(seq::GeneralSequence{A}) where {T,K,A<:DNAAlphabet}
    return make_kmer(Kmer{DNA,K}, seq)
end

function Kmer{T,K}(seq::GeneralSequence{A}) where {T,K,A<:RNAAlphabet}
    return make_kmer(Kmer{RNA,K}, seq)
end

Kmer{T,K}(x::Kmer{T,K}) where {T,K} = x
Kmer{T}(seq::AbstractString) where {T} = Kmer{T,length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A<:DNAAlphabet} = Kmer{DNA,length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A<:RNAAlphabet} = Kmer{RNA,length(seq)}(seq)
DNAKmer(seq::GeneralSequence{A}) where {A<:DNAAlphabet} = DNAKmer{length(seq)}(seq)
RNAKmer(seq::GeneralSequence{A}) where {A<:RNAAlphabet} = RNAKmer{length(seq)}(seq)

# create a kmer from a sequence whose elements are convertible to a nucleotide
function make_kmer(::Type{Kmer{T,K}}, seq) where {T,K}
    seqlen = length(seq)
    if seqlen > 32
        throw(ArgumentError("cannot create a k-mer loger than 32nt"))
    elseif seqlen != K
        throw(ArgumentError("cannot create a $(K)-mer from a sequence of length $(seqlen)"))
    end

    x = UInt64(0)
    for c in seq
        nt = convert(T, c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a k-mer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a k-mer with gaps"))
        end
        x = (x << 2) | UInt64(trailing_zeros(nt))
    end

    return Kmer{T,K}(x)
end

make_kmer(seq::NTuple{K,T}) where {K,T} = make_kmer(Kmer{T,K}, seq)

GeneralSequence(x::DNAKmer{K}) where {K} = DNASequence(x)
GeneralSequence(x::RNAKmer{K}) where {K} = RNASequence(x)
GeneralSequence{A}(x::DNAKmer{K}) where {A<:DNAAlphabet,K} = MutableBioSequence{A}([nt for nt in x])
GeneralSequence{A}(x::RNAKmer{K}) where {A<:RNAAlphabet,K} = MutableBioSequence{A}([nt for nt in x])
Base.convert(::Type{S}, seq::Kmer) where {S<:AbstractString} = S([Char(x) for x in seq])
Base.String(seq::Kmer) = convert(String, seq)


# Basic Functions
# ---------------

BioSymbols.alphabet(::Type{DNAKmer{k}}) where {k} = (DNA_A, DNA_C, DNA_G, DNA_T)
BioSymbols.alphabet(::Type{RNAKmer{k}}) where {k} = (RNA_A, RNA_C, RNA_G, RNA_U)
Alphabet(::Type{Kmer{T,N} where T<:NucleicAcid}) where N = Any
Alphabet(::Type{T}) where T <: DNAKmer = DNAAlphabet{2}
Alphabet(::Type{T}) where T <: RNAKmer = RNAAlphabet{2}

Base.hash(x::Kmer, h::UInt) = hash(UInt64(x), h)

kmersize(::Type{Kmer{T,k}}) where {T,k} = k
kmersize(kmer::Kmer) = kmersize(typeof(kmer))
Base.length(x::Kmer{T, K}) where {T,K} = kmersize(x)
Base.eltype(::Type{Kmer{T,k}}) where {T,k} = T

@inline function inbounds_getindex(x::Kmer{T,K}, i::Integer) where {T,K}
    return reinterpret(T, 0x01 << ((UInt64(x) >> 2(K - i)) & 0b11))
end

Base.summary(x::DNAKmer{k}) where {k} = string("DNA ", k, "-mer")
Base.summary(x::RNAKmer{k}) where {k} = string("RNA ", k, "-mer")

Base.:-(x::Kmer{T,K}, y::Integer) where {T,K} = Kmer{T,K}(UInt64(x) - y % UInt64)
Base.:+(x::Kmer{T,K}, y::Integer) where {T,K} = Kmer{T,K}(UInt64(x) + y % UInt64)
Base.:+(x::Integer, y::Kmer{T,K}) where {T,K} = y + x
Base.:(==)(x::Kmer{T,k}, y::Kmer{T,k}) where {T,k} = UInt64(x) == UInt64(y)
Base.isless(x::Kmer{T,K}, y::Kmer{T,K}) where {T,K} = isless(UInt64(x), UInt64(y))

function Base.typemin(::Type{Kmer{T,K}}) where {T,K}
    checkkmer(Kmer{T,K})
    return reinterpret(Kmer{T,K}, UInt64(0))
end

function Base.typemax(::Type{Kmer{T,K}}) where {T,K}
    checkkmer(Kmer{T,K})
    return reinterpret(Kmer{T,K}, ~UInt64(0) >> (64 - 2K))
end

@inline function checkkmer(::Type{Kmer{T,K}}) where {T,K}
    if !(1 ≤ K ≤ 32)
        throw(ArgumentError("the length K must be within 1..32"))
    end
end


# Other functions
# ---------------

"""
    complement(kmer::Kmer)

Return the complement of `kmer`.
"""
BioSymbols.complement(x::Kmer{T,k}) where {T,k} = Kmer{T,k}(~UInt64(x))

"""
    reverse(kmer::Kmer)

Return the reverse of `kmer`.
"""
Base.reverse(x::Kmer{T,k}) where {T,k} = Kmer{T,k}(nucrev2(UInt64(x)) >> (64 - 2k))

"""
    reverse_complement(kmer::Kmer)

Return the reverse complement of `kmer`
"""
reverse_complement(x::Kmer) = complement(reverse(x))

"""
    mismatches(a::Kmer, b::Kmer)

Return the number of mismatches between `a` and `b`.
"""
function mismatches(a::Kmer{T,k}, b::Kmer{T,k}) where {T,k}
    return count_nonzero_bitpairs(UInt64(a) ⊻ UInt64(b))
end

"""
    canonical(kmer::Kmer)

Return the canonical k-mer of `x`.

A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting k-mers in data that is not strand specific,
and thus observing k-mer is equivalent to observing its reverse complement.
"""
function canonical(x::Kmer)
    y = reverse_complement(x)
    return x < y ? x : y
end

function Base.rand(::Type{Kmer{T,k}}) where {T,k}
    return Kmer{T,k}(rand(UInt64))
end

function Base.rand(::Type{Kmer{T,k}}, size::Integer) where {T,k}
    return [rand(Kmer{T,k}) for _ in 1:size]
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


# Counters
# --------

function gc_content(kmer::Kmer{T,k}) where {T,k}
    if k == 0
        return 0.0
    else
        return (count_g(kmer) + count_c(kmer)) / k
    end
end

function count_a(kmer::Kmer{T,k}) where {T,k}
    return count_a(reinterpret(UInt64, kmer)) - (32 - k)
end

function count_c(kmer::Kmer{T,k}) where {T,k}
    return count_c(reinterpret(UInt64, kmer))
end

function count_g(kmer::Kmer{T,k}) where {T,k}
    return count_g(reinterpret(UInt64, kmer))
end

function count_t(kmer::Kmer{T,k}) where {T,k}
    return count_t(reinterpret(UInt64, kmer))
end

# Count A, C, T/U, G respectively in a kmer stored in a UInt64
function count_a(x::UInt64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::UInt64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::UInt64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::UInt64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)


# Shuffle
# -------

function Random.shuffle(kmer::Kmer{T,k}) where {T,k}
    # Fisher-Yates shuffle
    for i in 1:k-1
        j = rand(i:k)
        kmer = swap(kmer, i, j)
    end
    return kmer
end

# Swap two nucleotides at `i` and `j`.
function swap(kmer::Kmer{T,k}, i, j) where {T,k}
    i = 2k - 2i
    j = 2k - 2j
    b = convert(UInt64, kmer)
    x = ((b >> i) ⊻ (b >> j)) & UInt64(0x03)
    return Kmer{T,k}(b ⊻ ((x << i) | (x << j)))
end

# String literal
# --------------

macro kmer_str(seq)
    return DNAKmer(remove_newlines(seq))
end
