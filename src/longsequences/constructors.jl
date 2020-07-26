###
### Constructors
###
###
### Constructor methods for LongSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function LongSequence{A}(len::Integer) where {A<:Alphabet}
    if len < 0
        throw(ArgumentError("len must be non-negative"))
    end
    return LongSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), convert(Int, len))
end

LongSequence(::Type{DNA}) = LongDNASeq()
LongSequence(::Type{RNA}) = LongRNASeq()
LongSequence(::Type{AminoAcid}) = LongAminoAcidSeq()
LongSequence(::Type{Char}) = LongCharSeq()

function LongSequence()
    return LongSequence{VoidAlphabet}(Vector{UInt64}(), 0)
end

function LongSequence{A}(seq::LongSequence{A}, part::UnitRange) where A
    return seq[part]
end

function LongSequence{A}(s::Union{String, SubString{String}}) where {A<:Alphabet}
    return LongSequence{A}(s, codetype(A()))
end

# Generic method for String/Substring.
function LongSequence{A}(s::Union{String, SubString{String}}, ::AlphabetCode) where {A<:Alphabet}
    len = length(s)
    seq = LongSequence{A}(len)
    return copyto!(seq, 1, s, 1, len)
end

function LongSequence{A}(s::Union{String, SubString{String}}, ::AsciiAlphabet) where {A<:Alphabet}
    v = GC.@preserve s unsafe_wrap(Vector{UInt8}, pointer(s), ncodeunits(s))
    seq = LongSequence{A}(length(v))
    return encode_chunks!(seq, 1, v, 1, length(v))
end

function LongSequence{A}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = LongSequence{A}(len)
    return copyto!(seq, 1, src, startpos, len)
end

# create a subsequence
function LongSequence(other::LongSequence, part::UnitRange{<:Integer})
    checkbounds(other, part)
    subseq = typeof(other)(length(part))
    copyto!(subseq, 1, other, first(part), length(part))
    return subseq
end

function (::Type{T})(seq::BioSequence) where {T<:LongSequence}
    newseq = T(length(seq))
    @inbounds for i in eachindex(seq)
        newseq[i] = seq[i]
    end
    return newseq
end

function LongSequence(seq::BioSequence{A}) where {A <: Alphabet}
    return LongSequence{A}(seq)
end

function LongSequence{A}(seq::LongSequence{A}) where {A <: Alphabet}
    return LongSequence{A}(copy(seq.data), seq.len)
end

function (::Type{T})(seq::LongSequence{<:NucleicAcidAlphabet{N}}) where
         {N, T<:LongSequence{<:NucleicAcidAlphabet{N}}}
    return T(copy(seq.data), seq.len)
end

# This exists to fix ambiguity errors with Julia 1.0
function LongSequence{A}(seq::LongSequence{<:NucleicAcidAlphabet{N}}) where {N, A <: NucleicAcidAlphabet{N}}
	return LongSequence{A}(copy(seq.data), seq.len)
end

# Concatenate multiple sequences
function LongSequence{A}(chunks::LongSequence{A}...) where {A}
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = LongSequence{A}(len)
    offset = 1
    for chunk in chunks
        copyto!(seq, offset, chunk, 1, length(chunk))
        offset += length(chunk)
    end
    return seq
end

function Base.repeat(chunk::LongSequence{A}, n::Integer) where {A}
    seq = LongSequence{A}(length(chunk) * n)
    offset = 1
    for i in 1:n
        copyto!(seq, offset, chunk, 1, length(chunk))
        offset += length(chunk)
    end
    return seq
end

# Concatenation and Base.repeat operators
Base.:*(chunk::LongSequence{A}, chunks::LongSequence{A}...) where {A} =
    LongSequence{A}(chunk, chunks...)
Base.:^(chunk::LongSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::LongSequence{A}, len::Integer = length(seq)) where {A}
    return LongSequence{A}(len)
end
