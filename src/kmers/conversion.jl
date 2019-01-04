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

make_skipmer(seq::Tuple{DNA, Vararg{DNA, K}}) where {K} = make_skipmer(Skipmer{DNAAlphabet{2}, 2, 3, K + 1}, seq)
make_skipmer(seq::Tuple{RNA, Vararg{RNA, K}}) where {K} = make_skipmer(Skipmer{RNAAlphabet{2}, 2, 3, K + 1}, seq)
make_kmer(seq::Tuple{DNA, Vararg{DNA, K}}) where {K} = make_skipmer(Kmer{DNAAlphabet{2}, K + 1}, seq)
make_kmer(seq::Tuple{RNA, Vararg{RNA, K}}) where {K} = make_skipmer(Kmer{RNAAlphabet{2}, K + 1}, seq)
make_bigskipmer(seq::Tuple{DNA, Vararg{DNA, K}}) where {K} = make_skipmer(BigSkipmer{DNAAlphabet{2}, 2, 3, K + 1}, seq)
make_bigskipmer(seq::Tuple{RNA, Vararg{RNA, K}}) where {K} = make_skipmer(BigSkipmer{RNAAlphabet{2}, 2, 3, K + 1}, seq)
make_bigkmer(seq::Tuple{DNA, Vararg{DNA, K}}) where {K} = make_skipmer(BigKmer{DNAAlphabet{2}, K + 1}, seq)
make_bigkmer(seq::Tuple{RNA, Vararg{RNA, K}}) where {K} = make_skipmer(BigKmer{RNAAlphabet{2}, K + 1}, seq)

Skipmer(nts::T...) where {T<:NucleicAcid} = make_skipmer(nts)
Kmer(nts::T...) where {T<:NucleicAcid} = make_kmer(nts)
BigSkipmer(nts::T...) where {T<:NucleicAcid} = make_bigskipmer(nts)
BigKmer(nts::T...) where {T<:NucleicAcid} = make_bigkmer(nts)

function DNACodon(x::DNA, y::DNA, z::DNA)
    return make_skipmer(Kmer{DNAAlphabet{2}, 3}, (x, y, z))
end

function RNACodon(x::RNA, y::RNA, z::RNA)
    return make_skipmer(Kmer{RNAAlphabet{2},3}, (x, y, z))
end

@inline function make_mask(::Type{T}) where {T <: Union{Skipmer, BigSkipmer}}
    return (encoded_data_eltype(T)(1) << (2 * kmersize(T))) - 1
end

function Skipmer{A, M, N, K}(x::UInt64) where {A, M, N, K}
    checkskipmer(Skipmer{A, M, N, K})
    mask = make_mask(Skipmer{A, M, N, K})
    return reinterpret(Skipmer{A, M, N, K}, x & mask)
end
UInt64(x::Skipmer) = reinterpret(UInt64, x)
Base.convert(::Type{UInt64}, x::Skipmer) = reinterpret(UInt64, x)

function BigSkipmer{A, M, N, K}(x::UInt128) where {A, M, N, K}
    checkskipmer(BigSkipmer{A, M, N, K})
    mask = make_mask(BigSkipmer{A, M, N, K})
    return reinterpret(BigKmer{A, K}, x & mask)
end
UInt128(x::BigSkipmer) = reinterpret(UInt128, x)
Base.convert(::Type{UInt128}, x::BigSkipmer) = reinterpret(UInt128, x)

Skipmer{DNAAlphabet{2}, M, N, K}(x::Skipmer{RNAAlphabet{2}, M, N, K}) where {M, N, K} = reinterpret(Skipmer{DNAAlphabet{2}, M, N, K}, x)
Skipmer{RNAAlphabet{2}, M, N, K}(x::Skipmer{DNAAlphabet{2}, M, N, K}) where {M, N, K} = reinterpret(Skipmer{RNAAlphabet{2}, M, N, K}, x)
BigSkipmer{DNAAlphabet{2}, M, N, K}(x::BigSkipmer{RNAAlphabet{2}, M, N, K}) where {M, N, K} = reinterpret(BigSkipmer{DNAAlphabet{2}, M, N, K}, x)
BigSkipmer{RNAAlphabet{2}, M, N, K}(x::BigSkipmer{DNAAlphabet{2}, M, N, K}) where {M, N, K} = reinterpret(BigSkipmer{RNAAlphabet{2}, M, N, K}, x)

# Constructor for all concrete skiper types
function Skipmer{A, M, N, K}(seq::AbstractString) where {A, M, N, K}
    return make_skipmer(Skipmer{A, M, N, K}, seq)
end
Skipmer{A, M, N}(seq::AbstractString) where {A, M, N} = Skipmer{A, M, N, length(seq)}(seq)
Skipmer{A}(seq::AbstractString) where A = Skipmer{A, 2, 3, length(seq)}(seq)

function BigSkipmer{A, M, N, K}(seq::AbstractString) where {A, M, N, K}
    return make_skipmer(BigSkipmer{A, M, N, K}, seq)
end
BigSkipmer{A, M, N}(seq::AbstractString) where {A, M, N} = BigSkipmer{A, M, N, length(seq)}(seq)
BigSkipmer{A}(seq::AbstractString) where {A} = BigSkipmer{A, 2, 3, length(seq)}(seq)

# Constructor for all concrete skipmer types
function Skipmer{A, M, N, K}(seq::GeneralSequence) where {A, M, N, K}
    return make_skipmer(Skipmer{A, M, N, K}, seq)
end
Skipmer{A, M, N}(seq::GeneralSequence) where {A, M, N} = Skipmer{A, M, N, length(seq)}(seq)
Skipmer(seq::GeneralSequence{A}) where {A <: DNAAlphabet} = Skipmer{DNAAlphabet{2}, 2, 3, length(seq)}(seq)
Skipmer(seq::GeneralSequence{A}) where {A <: RNAAlphabet} = Skipmer{RNAAlphabet{2}, 2, 3, length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A <: DNAAlphabet} = Kmer{DNAAlphabet{2},length(seq)}(seq)
Kmer(seq::GeneralSequence{A}) where {A <: RNAAlphabet} = Kmer{RNAAlphabet{2},length(seq)}(seq)

function BigSkipmer{A, M, N, K}(seq::GeneralSequence) where {A, M, N, K}
    return make_skipmer(BigSkipmer{A, M, N, K}, seq)
end
BigSkipmer{A, M, N}(seq::GeneralSequence) where {A, M, N} = BigSkipmer{A, M, N, length(seq)}(seq)
BigSkipmer(seq::GeneralSequence{A}) where {A <: DNAAlphabet} = BigSkipmer{DNAAlphabet{2}, 2, 3, length(seq)}(seq)
BigSkipmer(seq::GeneralSequence{A}) where {A <: RNAAlphabet} = BigSkipmer{RNAAlphabet{2}, 2, 3, length(seq)}(seq)
BigKmer(seq::GeneralSequence{A}) where {A <: DNAAlphabet} = BigKmer{DNAAlphabet{2},length(seq)}(seq)
BigKmer(seq::GeneralSequence{A}) where {A <: RNAAlphabet} = BigKmer{RNAAlphabet{2},length(seq)}(seq)

Skipmer{A, M, N, K}(x::Skipmer{A, M, N, K}) where {A, M, N, K} = x

GeneralSequence(x::T) where {T <: Union{Skipmer{DNAAlphabet{2}}, BigSkipmer{DNAAlphabet{2}}}} = GeneralSequence{DNAAlphabet{2}}(x)
GeneralSequence(x::T) where {T <: Union{Skipmer{RNAAlphabet{2}}, BigSkipmer{RNAAlphabet{2}}}} = GeneralSequence{RNAAlphabet{2}}(x)
GeneralSequence{A}(x::T) where {A <: DNAAlphabet, T <: Union{Skipmer{DNAAlphabet{2}}, BigSkipmer{DNAAlphabet{2}}}} = GeneralSequence{A}([nt for nt in x])
GeneralSequence{A}(x::T) where {A <: RNAAlphabet, T <: Union{Skipmer{RNAAlphabet{2}}, BigSkipmer{RNAAlphabet{2}}}} = GeneralSequence{A}([nt for nt in x])