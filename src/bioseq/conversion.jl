# Conversion & Promotion
# ======================
#
# Conversion methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Promotion
# ---------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{BioSequence{A}}, ::Type{BioSequence{B}}) where {A<:$alph,B<:$alph}
        return BioSequence{promote_rule(A,B)}
    end
end

# Conversion
# ----------

# Conversion between sequences of different alphabet size.

function _generic_convert(::Type{T}, seq) where T <: Sequence
    # TODO: make it faster with bit-parallel algorithm
    newseq = T(length(seq))
    for (i, x) in enumerate(seq)
        unsafe_setindex!(newseq, x, i)
    end
    return newseq
end

function convert(::Type{BioSequence{DNAAlphabet{2}}}, seq::BioSequence{DNAAlphabet{4}})
    return _generic_convert(DNAAlphabet{2}, seq)
end
BioSequence{DNAAlphabet{2}}(seq::BioSequence{DNAAlphabet{4}}) = convert(BioSequence{DNAAlphabet{2}}, seq)

function convert(::Type{BioSequence{DNAAlphabet{4}}}, seq::BioSequence{DNAAlphabet{2}})
    return _generic_convert(DNAAlphabet{4}, seq)
end
BioSequence{DNAAlphabet{4}}(seq::BioSequence{DNAAlphabet{2}}) = convert(BioSequence{DNAAlphabet{4}}, seq)

# Conversion between DNA and RNA sequences.
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)], n in (2, 4)
    # NOTE: assumes that binary representation is identical between DNA and RNA
    @eval function BioSequence{$(A1{n})}(seq::BioSequence{$(A2{n})})
        newseq = BioSequence{$(A1{n})}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
    @eval function Base.convert(::Type{BioSequence{$(A1{n})}}, seq::BioSequence{$(A2{n})})
        return BioSequence{$(A1{n})}(seq)
    end
end

# Convert from a DNA or RNA vector to a BioSequence.
function BioSequence{A}(seq::AbstractVector{DNA}) where {A<:DNAAlphabet}
    return BioSequence{A}(seq, 1, lastindex(seq))
end
function BioSequence{A}(seq::AbstractVector{RNA}) where {A<:RNAAlphabet}
    return BioSequence{A}(seq, 1, lastindex(seq))
end
function AminoAcidSequence(seq::AbstractVector{AminoAcid})
    return AminoAcidSequence(seq, 1, lastindex(seq))
end

# Convert from a BioSequence to to a DNA or RNA vector
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
function Base.convert(::Type{Vector{DNA}}, seq::BioSequence{<:DNAAlphabet})
    return collect(seq)
end
function Base.convert(::Type{Vector{RNA}}, seq::BioSequence{<:RNAAlphabet})
    return collect(seq)
end
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

# Covert from a string to a BioSequence and _vice versa_.
function Base.convert(::Type{S}, seq::BioSequence) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end
Base.String(seq::BioSequence) = convert(String, seq)
Base.convert(::Type{BioSequence{A}}, seq::S) where {S<:AbstractString,A} = BioSequence{A}(seq)

