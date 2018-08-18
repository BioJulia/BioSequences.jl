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
        return BioSequence{promote_rule(A, B)}
    end
end

# Conversion
# ----------

function Base.convert(::Type{BioSequence{DNAAlphabet{4}}}, seq::BioSequence{DNAAlphabet{2}})
    return BioSequence{DNAAlphabet{4}}(seq)
end

function Base.convert(::Type{BioSequence{DNAAlphabet{2}}}, seq::BioSequence{DNAAlphabet{4}})
    return BioSequence{DNAAlphabet{2}}(seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{4}}}, seq::BioSequence{RNAAlphabet{2}})
    return BioSequence{RNAAlphabet{4}}(seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{2}}}, seq::BioSequence{RNAAlphabet{4}})
    return BioSequence{RNAAlphabet{2}}(seq)
end

function Base.convert(::Type{BioSequence{DNAAlphabet{2}}}, seq::BioSequence{RNAAlphabet{2}})
    return BioSequence{DNAAlphabet{2}}(seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{2}}}, seq::BioSequence{DNAAlphabet{2}})
    return BioSequence{RNAAlphabet{2}}(seq)
end

function Base.convert(::Type{BioSequence{DNAAlphabet{4}}}, seq::BioSequence{RNAAlphabet{4}})
    return BioSequence{DNAAlphabet{4}}(seq)
end

function Base.convert(::Type{BioSequence{RNAAlphabet{4}}}, seq::BioSequence{DNAAlphabet{4}})
    return BioSequence{RNAAlphabet{4}}(seq)
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
Base.convert(::Type{BioSequence{A}}, seq::AbstractString) where A = BioSequence{A}(seq)

