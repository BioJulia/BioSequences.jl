###
### Conversion & Promotion
###
###
### Conversion methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Create a lookup table from biosymbol to the UInt8 for the character that would
# represent it in a string, e.g. DNA_G -> UInt8('G')
for alphabettype in ("DNA", "RNA", "AminoAcid")
    tablename = Symbol(uppercase(alphabettype), "_TO_BYTE")
    typ = Symbol(alphabettype)
    @eval begin
        const $(tablename) = let
            alph = alphabet($(typ))
            bytes = zeros(UInt8, length(alph))
            @inbounds for letter in alph
                bytes[reinterpret(UInt8, letter) + 1] = UInt8(Char(letter))
            end
            Tuple(bytes)
        end
        stringbyte(x::$(typ)) = @inbounds $(tablename)[reinterpret(UInt8, x) + 1]
    end
end

# Less efficient fallback. Should only be called for symbols of AsciiAlphabet
stringbyte(x::BioSymbol) = UInt8(Char(x))

abstract type AlphabetCode end
struct AsciiAlphabet <: AlphabetCode end
struct UnicodeAlphabet <: AlphabetCode end

function codetype(::A) where {A <: Union{DNAAlphabet{2}, DNAAlphabet{4},
                                         RNAAlphabet{2}, RNAAlphabet{4},
                                         AminoAcidAlphabet}}
    return AsciiAlphabet()
end
codetype(::Alphabet) = UnicodeAlphabet()

function Base.convert(::Type{S}, seq::BioSequence) where {S<:AbstractString}
    return convert(S, seq, codetype(Alphabet(seq)))
end

@inline function Base.convert(::Type{S}, seq::BioSequence, ::AlphabetCode) where {S<:AbstractString}
    return S([Char(x) for x in seq])
end

@inline function Base.convert(::Type{String}, seq::BioSequence, ::AsciiAlphabet)
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
