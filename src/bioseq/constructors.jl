# Constructors
# ============
#
# Constructor methods for Biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function BioSequence{A}(len::Integer) where {A<:Alphabet}
    return BioSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), 1:len, false)
end

BioSequence(::Type{DNA}) = DNASequence()
BioSequence(::Type{RNA}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()
BioSequence(::Type{Char}) = CharSequence()

function BioSequence()
    return BioSequence{VoidAlphabet}(Vector{UInt64}(), 0:-1, false)
end

function BioSequence{A}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = BioSequence{A}(len)
    return encode_copy!(seq, 1, src, startpos, len)
end

# create a subsequence
function BioSequence(other::BioSequence{A}, part::UnitRange{<:Integer}) where {A}
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end

function BioSequence{A}(other::BioSequence{A}, part::UnitRange) where {A}
    return BioSequence(other, part)
end

# concatenate chunks
function BioSequence{A}(chunks::BioSequence{A}...) where {A}
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = BioSequence{A}(len)
    offset = 1
    for chunk in chunks
        copyto!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

function Base.repeat(chunk::BioSequence{A}, n::Integer) where {A}
    seq = BioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copyto!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# operators for concat and repeat
Base.:*(chunk::BioSequence{A}, chunks::BioSequence{A}...) where {A} =
    BioSequence{A}(chunk, chunks...)
Base.:^(chunk::BioSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::BioSequence{A}, len::Integer=length(seq)) where {A}
    return BioSequence{A}(len)
end
