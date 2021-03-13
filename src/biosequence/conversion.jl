###
### Conversion & Promotion
###
###
### Conversion methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.convert(::Type{S}, seq::BioSequence) where {S<:AbstractString}
    return convert(S, seq, codetype(Alphabet(seq)))
end

@inline function Base.convert(::Type{S}, seq::BioSequence, ::AlphabetCode) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end

function Base.convert(::Type{String}, seq::BioSequence, ::AsciiAlphabet)
    len = length(seq)
    str = Base._string_n(len)
    GC.@preserve str begin
        p = pointer(str)
        @inbounds for i in 1:len
            unsafe_store!(p, stringbyte(seq[i]), i)
        end
    end
    return str
end

Base.String(seq::BioSequence) = convert(String, seq)
Base.convert(::Type{Vector{DNA}}, seq::BioSequence{<:DNAAlphabet}) = collect(seq)
Base.convert(::Type{Vector{RNA}}, seq::BioSequence{<:RNAAlphabet}) = collect(seq)
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSeq) = collect(seq)
