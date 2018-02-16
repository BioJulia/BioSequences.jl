# Constructors
# ============
#
# Constructor methods for Biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function MutableBioSequence{A}(len::Integer) where {A<:Alphabet}
    return MutableBioSequence{A}(Vector{UInt64}(seq_data_len(A, len)), 1:len, false)
end

MutableBioSequence(::Type{DNA}) = DNASequence()
MutableBioSequence(::Type{RNA}) = RNASequence()
MutableBioSequence(::Type{AminoAcid}) = AminoAcidSequence()
MutableBioSequence(::Type{Char}) = CharSequence()

function MutableBioSequence()
    return MutableBioSequence{VoidAlphabet}(Vector{UInt64}(0), 0:-1, false)
end

function MutableBioSequence{A}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = MutableBioSequence{A}(len)
    #println("Made empty sequence ", seq)
    #println("Making the encode_copy!")
    return encode_copy!(seq, 1, src, startpos, len)
end

# create a subsequence
function MutableBioSequence(other::MutableBioSequence{A}, part::UnitRange{<:Integer}) where {A}
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = MutableBioSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end

function MutableBioSequence{A}(other::MutableBioSequence{A}, part::UnitRange) where {A}
    return MutableBioSequence(other, part)
end

# concatenate chunks
function MutableBioSequence{A}(chunks::MutableBioSequence{A}...) where {A}
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = MutableBioSequence{A}(len)
    offset = 1
    for chunk in chunks
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

function Base.repeat(chunk::MutableBioSequence{A}, n::Integer) where {A}
    seq = MutableBioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# operators for concat and repeat
Base.:*(chunk::MutableBioSequence{A}, chunks::MutableBioSequence{A}...) where {A} =
    MutableBioSequence{A}(chunk, chunks...)
Base.:^(chunk::MutableBioSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::MutableBioSequence{A}, len::Integer=length(seq)) where {A}
    return MutableBioSequence{A}(len)
end
