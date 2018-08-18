# Constructors
# ============
#
# Constructor methods for Biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Construct a blank sequence of length `len`.
function BioSequence{A}(len::Integer) where {A<:Alphabet}
    return BioSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), 1:len, false)
end

# Construct a void sequence.
function BioSequence()
    return BioSequence{VoidAlphabet}(Vector{UInt64}(), 0:-1, false)
end

BioSequence(::Type{DNA}) = DNASequence()
BioSequence(::Type{RNA}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()
BioSequence(::Type{Char}) = CharSequence()

# Create a sequence from a string or a vector.
function BioSequence{A}(
        src::Union{AbstractString, AbstractVector},
        startpos::Integer = 1,
        stoppos::Integer = length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = BioSequence{A}(len)
    return encode_copy!(seq, 1, src, startpos, len)
end

# Create a sequence that is a subsequence of another sequence. 
function BioSequence{A}(
        other::BioSequence{A},
        part::UnitRange{<:Integer}) where {A<:Alphabet}
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end
function BioSequence(other::BioSequence{A}, part::UnitRange) where {A<:Alphabet}
    return BioSequence{A}(other, part)
end

# Concatenate multiple sequences.
function BioSequence{A}(chunks::BioSequence{A}...) where {A<:Alphabet}
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

# Create a 4 bit DNA/RNA sequence from a 2 bit DNA/RNA sequence, and vice-versa.
function BioSequence{DNAAlphabet{4}}(seq::BioSequence{DNAAlphabet{2}})
    newseq = BioSequence{DNAAlphabet{4}}(length(seq))
    for (i, x) in enumerate(seq)
        unsafe_setindex!(newseq, x, i)
    end
    return newseq
end
function BioSequence{DNAAlphabet{2}}(seq::BioSequence{DNAAlphabet{4}})
    newseq = BioSequence{DNAAlphabet{2}}(length(seq))
    for (i, x) in enumerate(seq)
        unsafe_setindex!(newseq, x, i)
    end
    return newseq
end
function BioSequence{RNAAlphabet{4}}(seq::BioSequence{RNAAlphabet{2}})
    newseq = BioSequence{RNAAlphabet{4}}(length(seq))
    for (i, x) in enumerate(seq)
        unsafe_setindex!(newseq, x, i)
    end
    return newseq
end
function BioSequence{RNAAlphabet{2}}(seq::BioSequence{RNAAlphabet{4}})
    newseq = BioSequence{RNAAlphabet{2}}(length(seq))
    for (i, x) in enumerate(seq)
        unsafe_setindex!(newseq, x, i)
    end
    return newseq
end

function Base.repeat(chunk::BioSequence{A}, n::Integer) where {A<:Alphabet}
    seq = BioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copyto!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# Concatenation and Base.repeat operators. 
Base.:*(chunk::BioSequence{A}, chunks::BioSequence{A}...) where {A} =
    BioSequence{A}(chunk, chunks...)
Base.:^(chunk::BioSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::BioSequence{A}, len::Integer = length(seq)) where {A<:Alphabet}
    return BioSequence{A}(len)
end


