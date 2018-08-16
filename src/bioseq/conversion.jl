# Conversion & Promotion
# ======================
#
# Conversion methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Promotion
# ---------

function Base.promote_rule(::Type{BioSequence{A}}, ::Type{BioSequence{B}}) where {A<:Alphabet,B<:Alphabet}
    return BioSequence{promote_rule(A, B)}
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

for A in (DNAAlphabet, RNAAlphabet)
    @eval function Base.convert(::Type{BioSequence{$A{2}}}, seq::BioSequence{$A{4}})
        return _generic_convert(BioSequence{$A{2}}, seq)
    end
    @eval function Base.convert(::Type{BioSequence{$A{4}}}, seq::BioSequence{$A{2}})
        return _generic_convert(BioSequence{$A{4}}, seq)
    end
end

#=
function BioSequence{DNAAlphabet{2}}(seq::BioSequence{DNAAlphabet{4}})
    return convert(BioSequence{DNAAlphabet{2}}, seq)
end
=#


#BioSequence{DNAAlphabet{4}}(seq::BioSequence{DNAAlphabet{2}}) = convert(BioSequence{DNAAlphabet{4}}, seq)

# Conversion between DNA and RNA sequences.
function _same_bits_convert(::Type{T}, seq) where T <: Sequence
    newseq = T(seq.data, seq.part, true)
    seq.shared = true
    return newseq
end

function Base.convert(::Type{BioSequence{DNAAlphabet{2}}}, seq::BioSequence{RNAAlphabet{2}})
    return _same_bits_convert(BioSequence{DNAAlphabet{2}}, seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{2}}}, seq::BioSequence{DNAAlphabet{2}})
    return _same_bits_convert(BioSequence{RNAAlphabet{2}}, seq)
end

function Base.convert(::Type{BioSequence{DNAAlphabet{4}}}, seq::BioSequence{RNAAlphabet{4}})
    return _same_bits_convert(BioSequence{DNAAlphabet{4}}, seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{4}}}, seq::BioSequence{DNAAlphabet{4}})
    return _same_bits_convert(BioSequence{RNAAlphabet{4}}, seq)
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

