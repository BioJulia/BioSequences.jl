###
### Conversion & Promotion
###
###
### Conversion methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.convert(::Type{S}, seq::BioSequence) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end
Base.String(seq::BioSequence) = convert(String, seq)

Base.convert(::Type{Vector{DNA}}, seq::BioSequence{<:DNAAlphabet}) = collect(seq)
Base.convert(::Type{Vector{RNA}}, seq::BioSequence{<:RNAAlphabet}) = collect(seq)
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSeq) = collect(seq)