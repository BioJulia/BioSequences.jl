###
### Conversion & Promotion
###
###
### Conversion methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function (::Type{S})(seq::BioSequence) where {S <: AbstractString}
    _string(S, seq, codetype(Alphabet(seq)))
end

function _string(::Type{S}, seq::BioSequence, ::AlphabetCode) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end

function _string(::Type{String}, seq::BioSequence, ::AsciiAlphabet)
    String([stringbyte(s) for s in seq])
end
