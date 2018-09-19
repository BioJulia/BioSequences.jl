# Conversion & Construction
# -------------------------

# Create a skipmer from a sequence whose elements are convertible to a nucleotide
function make_skipmer(::Type{T}, seq) where {T <: Union{Skipmer, BigSkipmer}}
    seqlen = length(seq)
    if seqlen != kmersize(T)
        throw(ArgumentError("seq does not contain the correct number of nucleotides ($seqlen ≠ $(kmersize(T)))"))
    end

    x = encoded_data_eltype(T)(0)
    for c in seq
        nt = convert(eltype(T), c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a skipmer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a skipmer with gaps"))
        end
        x = (x << 2) | encoded_data_eltype(T)(twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end

    return T(x)
end

make_skipmer(seq::NTuple{K,T}) where {K,T} = make_skipmer(Kmer{T,K}, seq)

function Kmer(nts::T...) where {T<:NucleicAcid}
    return make_skipmer(nts)
end

function DNACodon(x::DNA, y::DNA, z::DNA)
    return make_skipmer((x, y, z))
end

function RNACodon(x::RNA, y::RNA, z::RNA)
    return make_skipmer((x, y, z))
end

@inline function make_mask(::Type{T}) where {T <: Union{Skipmer, BigSkipmer}}
    return (encoded_data_eltype(T)(1) << (2 * kmersize(T))) - 1
end

function Skipmer{T, M, N, K}(x::UInt64) where {T, M, N, K}
    checkskipmer(Skipmer{T, M, N, K})
    mask = make_mask(Skipmer{T, M, N, K})
    return reinterpret(Skipmer{T, M, N, K}, x & mask)
end
UInt64(x::Skipmer) = reinterpret(UInt64, x)
Base.convert(::Type{UInt64}, x::Skipmer) = reinterpret(UInt64, x)

function BigSkipmer{T, M, N, K}(x::UInt128) where {T, M, N, K}
    checkskipkmer(BigSkipmer{T, M, N, K})
    mask = make_mask(BigSkipmer{T, M, N, K})
    return reinterpret(BigKmer{T, K}, x & mask)
end
UInt128(x::BigSkipmer) = reinterpret(UInt128, x)
Base.convert(::Type{UInt128}, x::BigSkipmer) = reinterpret(UInt128, x)

Skipmer{DNA, M, N, K}(x::Skipmer{RNA, M, N, K}) where {M, N, K} = reinterpret(Skipmer{DNA, M, N, K}, x)
Skipmer{RNA, M, N, K}(x::Skipmer{DNA, M, N, K}) where {M, N, K} = reinterpret(Skipmer{RNA, M, N, K}, x)
BigSkipmer{DNA, M, N, K}(x::BigSkipmer{RNA, M, N, K}) where {M, N, K} = reinterpret(BigSkipmer{DNA, M, N, K}, x)
BigSkipmer{RNA, M, N, K}(x::BigSkipmer{DNA, M, N, K}) where {M, N, K} = reinterpret(BigSkipmer{RNA, M, N, K}, x)

# Constructor for all concrete skiper types
function Skipmer{T, M, N, K}(seq::AbstractString) where {T, M, N, K}
    return make_skipmer(Skipmer{T, M, N, K}, seq)
end
Skipmer{T, M, N}(seq::AbstractString) where {T, M, N} = Skipmer{T, M, N, length(seq)}(seq)
Skipmer{T}(seq::AbstractString) where T = Skipmer{T, 2, 3, length(seq)}(seq)

# Constructor for all concrete skipmer types
function Skipmer{T, M, N, K}(seq::GeneralSequence) where {T, M, N, K}
    return BioSequences.make_skipmer(Skipmer{T, M, N, K}, seq)
end
Skipmer{T, M, N}(seq::GeneralSequence) where {T, M, N} = Skipmer{T, M, N, length(seq)}(seq)
Skipmer(seq::GeneralSequence) = Skipmer{eltype(seq), 2, 3, length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A <: NucleicAcidAlphabet} = Kmer{eltype(A),length(seq)}(seq)

Skipmer{T, M, N, K}(x::Skipmer{T, M, N, K}) where {T, M, N, K} = x

GeneralSequence(x::Skipmer{DNA, M, N, K}) where {M, N, K} = GeneralSequence{DNAAlphabet{2}}(x)
GeneralSequence(x::Skipmer{RNA, M, N, K}) where {M, N, K} = GeneralSequence{RNAAlphabet{2}}(x)
GeneralSequence{A}(x::Skipmer{DNA, M, N, K}) where {A <: DNAAlphabet, M, N, K} = GeneralSequence{A}([nt for nt in x])
GeneralSequence{A}(x::Skipmer{RNA, M, N, K}) where {A <: RNAAlphabet, M, N, K} = GeneralSequence{A}([nt for nt in x])