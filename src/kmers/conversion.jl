# Conversion & Construction
# -------------------------

# Create a skipmer from a sequence whose elements are convertible to a nucleotide
function make_skipmer(::Type{T}, seq) where {T <: Union{Skipmer, BigSkipmer}}
    checkskipmer(T)
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

make_skipmer(seq::Tuple{T, Vararg{T, K}}) where {K, T<:NucleicAcid} = make_skipmer(Skipmer{T, 2, 3, K + 1}, seq)
make_kmer(seq::Tuple{T, Vararg{T, K}}) where {K, T<:NucleicAcid} = make_skipmer(Kmer{T, K + 1}, seq)
make_bigskipmer(seq::Tuple{T, Vararg{T, K}}) where {K, T<:NucleicAcid} = make_skipmer(BigSkipmer{T, 2, 3, K + 1}, seq)
make_bigkmer(seq::Tuple{T, Vararg{T, K}}) where {K, T<:NucleicAcid} = make_skipmer(BigKmer{T, K + 1}, seq)

Skipmer(nts::T...) where {T<:NucleicAcid} = make_skipmer(nts)
Kmer(nts::T...) where {T<:NucleicAcid} = make_kmer(nts)
BigSkipmer(nts::T...) where {T<:NucleicAcid} = make_bigskipmer(nts)
BigKmer(nts::T...) where {T<:NucleicAcid} = make_bigkmer(nts)

function DNACodon(x::DNA, y::DNA, z::DNA)
    return make_skipmer(Kmer{DNA, 3}, (x, y, z))
end

function RNACodon(x::RNA, y::RNA, z::RNA)
    return make_skipmer(Kmer{RNA,3}, (x, y, z))
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
    checkskipmer(BigSkipmer{T, M, N, K})
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

function BigSkipmer{T, M, N, K}(seq::AbstractString) where {T, M, N, K}
    return make_skipmer(BigSkipmer{T, M, N, K}, seq)
end
BigSkipmer{T, M, N}(seq::AbstractString) where {T, M, N} = BigSkipmer{T, M, N, length(seq)}(seq)
BigSkipmer{T}(seq::AbstractString) where {T} = BigSkipmer{T, 2, 3, length(seq)}(seq)

# Constructor for all concrete skipmer types
function Skipmer{T, M, N, K}(seq::GeneralSequence) where {T, M, N, K}
    return make_skipmer(Skipmer{T, M, N, K}, seq)
end
Skipmer{T, M, N}(seq::GeneralSequence) where {T, M, N} = Skipmer{T, M, N, length(seq)}(seq)
Skipmer(seq::GeneralSequence) = Skipmer{eltype(seq), 2, 3, length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A <: NucleicAcidAlphabet} = Kmer{eltype(A),length(seq)}(seq)

function BigSkipmer{T, M, N, K}(seq::GeneralSequence) where {T, M, N, K}
    return make_skipmer(BigSkipmer{T, M, N, K}, seq)
end
BigSkipmer{T, M, N}(seq::GeneralSequence) where {T, M, N} = BigSkipmer{T, M, N, length(seq)}(seq)
BigSkipmer(seq::GeneralSequence) = BigSkipmer{eltype(seq), 2, 3, length(seq)}(seq)
BigKmer(seq::GeneralSequence{A}) where {A <: NucleicAcidAlphabet} = BigKmer{eltype(A),length(seq)}(seq)

Skipmer{T, M, N, K}(x::Skipmer{T, M, N, K}) where {T, M, N, K} = x

GeneralSequence(x::T) where {T <: Union{Skipmer{DNA}, BigSkipmer{DNA}}} = GeneralSequence{DNAAlphabet{2}}(x)
GeneralSequence(x::T) where {T <: Union{Skipmer{RNA}, BigSkipmer{RNA}}} = GeneralSequence{RNAAlphabet{2}}(x)
GeneralSequence{A}(x::T) where {A <: DNAAlphabet, T <: Union{Skipmer{DNA}, BigSkipmer{DNA}}} = GeneralSequence{A}([nt for nt in x])
GeneralSequence{A}(x::T) where {A <: RNAAlphabet, T <: Union{Skipmer{RNA}, BigSkipmer{RNA}}} = GeneralSequence{A}([nt for nt in x])