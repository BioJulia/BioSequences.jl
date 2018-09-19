# Constructors
# ============
#
# Constructor methods for Biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function GeneralSequence{A}(len::Integer) where {A<:Alphabet}
    return GeneralSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), 1:len, false)
end

GeneralSequence(::Type{DNA}) = DNASequence()
GeneralSequence(::Type{RNA}) = RNASequence()
GeneralSequence(::Type{AminoAcid}) = AminoAcidSequence()
GeneralSequence(::Type{Char}) = CharSequence()

function GeneralSequence()
    return GeneralSequence{VoidAlphabet}(Vector{UInt64}(0), 0:-1, false)
end

function GeneralSequence{A}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = GeneralSequence{A}(len)
    #println("Made empty sequence ", seq)
    #println("Making the encode_copy!")
    return encode_copy!(seq, 1, src, startpos, len)
end

# create a subsequence
function GeneralSequence(other::GeneralSequence{A}, part::UnitRange{<:Integer}) where {A}
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = GeneralSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end

function GeneralSequence{A}(other::GeneralSequence{A}, part::UnitRange) where {A}
    return GeneralSequence(other, part)
end

# Create a 4 bit DNA/RNA sequences from a 2 bit DNA/RNA sequences, and vice-versa.
for (alpha, alphb) in [(DNAAlphabet{4}, DNAAlphabet{2}), # DNA to DNA
                       (DNAAlphabet{2}, DNAAlphabet{4}),
                       (RNAAlphabet{4}, RNAAlphabet{2}), # RNA to RNA
                       (RNAAlphabet{2}, RNAAlphabet{4}),
                       (DNAAlphabet{2}, RNAAlphabet{4}), # DNA to RNA
                       (DNAAlphabet{4}, RNAAlphabet{2}),
                       (RNAAlphabet{4}, DNAAlphabet{2}), # RNA to DNA
                       (RNAAlphabet{2}, DNAAlphabet{4})]
    
    @eval function (::Type{GeneralSequence{$alpha}})(seq::GeneralSequence{$alphb})
        newseq = GeneralSequence{$alpha}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end
end

for (alpha, alphb) in [(DNAAlphabet{2}, RNAAlphabet{2}),
                       (RNAAlphabet{2}, DNAAlphabet{2}),
                       (DNAAlphabet{4}, RNAAlphabet{4}),
                       (RNAAlphabet{4}, DNAAlphabet{4})]
    
    @eval function (::Type{GeneralSequence{$alpha}})(seq::GeneralSequence{$alphb})
        newseq = GeneralSequence{$alpha}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
end


# Concatenate multiple sequences
function GeneralSequence{A}(chunks::GeneralSequence{A}...) where {A}
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = GeneralSequence{A}(len)
    offset = 1
    for chunk in chunks
        copyto!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

function Base.repeat(chunk::GeneralSequence{A}, n::Integer) where {A}
    seq = GeneralSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copyto!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# Concatenation and Base.repeat operators
Base.:*(chunk::GeneralSequence{A}, chunks::GeneralSequence{A}...) where {A} =
    GeneralSequence{A}(chunk, chunks...)
Base.:^(chunk::GeneralSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::GeneralSequence{A}, len::Integer = length(seq)) where {A}
    return GeneralSequence{A}(len)
end
