# Conversion & Promotion
# ======================
#
# Conversion methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.convert(::Type{S}, seq::BioSequence) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end
Base.String(seq::BioSequence) = convert(String, seq)