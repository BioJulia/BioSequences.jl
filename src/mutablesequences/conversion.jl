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
    @eval function Base.promote_rule(::Type{GeneralSequence{A}}, ::Type{GeneralSequence{B}}) where {A<:$alph,B<:$alph}
        return GeneralSequence{promote_rule(A, B)}
    end
end

# Conversion
# ----------

# Create a 4 bit DNA/RNA sequence from a 2 bit DNA/RNA sequence, and vice-versa.
for (alpha, alphb) in [(DNAAlphabet{4}, DNAAlphabet{2}), # DNA to DNA
                       (DNAAlphabet{2}, DNAAlphabet{4}),
                       (RNAAlphabet{4}, RNAAlphabet{2}), # RNA to RNA
                       (RNAAlphabet{2}, RNAAlphabet{4}),
                       (DNAAlphabet{2}, RNAAlphabet{4}), # DNA to RNA
                       (DNAAlphabet{4}, RNAAlphabet{2}),
                       (DNAAlphabet{2}, RNAAlphabet{2}),
                       (DNAAlphabet{4}, RNAAlphabet{4}),
                       (RNAAlphabet{4}, DNAAlphabet{2}), # RNA to DNA
                       (RNAAlphabet{2}, DNAAlphabet{4}),
                       (RNAAlphabet{2}, DNAAlphabet{2}),
                       (RNAAlphabet{4}, DNAAlphabet{4})]
    
    @eval function Base.convert(::Type{GeneralSequence{$alpha}}, seq::GeneralSequence{$alphb})
        return GeneralSequence{$alpha}(seq)
    end
end

# Convert from a GeneralSequence to to a DNA or RNA vector
function Base.convert(::Type{GeneralSequence{A}}, seq::Vector) where A<:Alphabet
    return GeneralSequence{A}(seq)
end
Base.convert(::Type{Vector}, seq::GeneralSequence) = collect(seq)
function Base.convert(::Type{Vector{DNA}}, seq::GeneralSequence{<:DNAAlphabet})
    return collect(seq)
end
function Base.convert(::Type{Vector{RNA}}, seq::GeneralSequence{<:RNAAlphabet})
    return collect(seq)
end
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

#= Covert from a string to a BioSequence and _vice versa_.
function Base.convert(::Type{S}, seq::GeneralSequence) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end
Base.String(seq::GeneralSequence) = convert(String, seq)
=#
Base.convert(::Type{GeneralSequence{A}}, seq::AbstractString) where A = GeneralSequence{A}(seq)

