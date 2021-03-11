###
### Conversion & Construction
###
###
### Conversion methods for Mers.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Create a Mer from a sequence whose elements are convertible to a nucleotide.
function (::Type{T})(seq) where {T<:AbstractMer}
    checkmer(T)
    seqlen = length(seq)
    if seqlen != ksize(T)
        throw(ArgumentError("seq does not contain the correct number of nucleotides ($seqlen ≠ $(ksize(T)))"))
    end
    if seqlen > capacity(T)
        throw(ArgumentError("Cannot build a mer longer than $(capacity(T))bp long"))
    end

    x = zero(encoded_data_type(T))
    for c in seq
        nt = convert(eltype(T), c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a mer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a mer with gaps"))
        end
        x = (x << 2) | encoded_data_type(T)(twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end

    return reinterpret(T, x)
end

(::Type{Mer{A}})(seq) where {A} = Mer{A,length(seq)}(seq)
(::Type{BigMer{A}})(seq) where {A} = BigMer{A,length(seq)}(seq)

function (::Type{T})(nts::Vararg{DNA,K}) where {K,T<:AbstractMer{DNAAlphabet{2},K}}
    return T(nts)
end
DNAMer(nts::Vararg{DNA,K}) where {K} = DNAMer{K}(nts)
BigDNAMer(nts::Vararg{DNA,K}) where {K} = BigDNAMer{K}(nts)

function (::Type{T})(nts::Vararg{RNA,K}) where {K,T<:AbstractMer{RNAAlphabet{2},K}}
    return T(nts)
end
RNAMer(nts::Vararg{RNA,K}) where {K} = RNAMer{K}(nts)
BigRNAMer(nts::Vararg{RNA,K}) where {K} = BigRNAMer{K}(nts)

LongSequence{A}(x::AbstractMer{DNAAlphabet{2},K}) where {A<:DNAAlphabet,K} = LongSequence{A}([nt for nt in x])
LongSequence{A}(x::AbstractMer{RNAAlphabet{2},K}) where {A<:RNAAlphabet,K} = LongSequence{A}([nt for nt in x])
LongSequence(x::AbstractMer{A,K}) where {A,K} = LongSequence{A}([nt for nt in x])
