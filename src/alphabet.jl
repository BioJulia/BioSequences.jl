# Alphabet
# ========
#
# Alphabet of biological symbols.
#
# Subtypes of `Alphabet` represent a domain of biological characters. For
# example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides (i.e. A, C,
# G, and T). These types are used for parameterizing biological sequences and so
# on. A pair of encoder and decoder is associated with an alphabet, which maps
# values between binary and Julia-level representation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Alphabet of biological characters.
"""
abstract type Alphabet end

"""
DNA nucleotide alphabet.
"""
struct DNAAlphabet{n} <: Alphabet end

"""
RNA nucleotide alphabet.
"""
struct RNAAlphabet{n} <: Alphabet end

"""
Amino acid alphabet.
"""
struct AminoAcidAlphabet <: Alphabet end

"""
General character alphabet.
"""
struct CharAlphabet <: Alphabet end

"""
Void alphabet (internal use only).
"""
struct VoidAlphabet <: Alphabet end

const NucAlphs = Union{DNAAlphabet,RNAAlphabet}
const TwoBitNucs = Union{DNAAlphabet{2}, RNAAlphabet{2}}
const FourBitNucs = Union{DNAAlphabet{4}, RNAAlphabet{4}}

"""
The number of bits to represent the alphabet.
"""
function bitsof end

for n in (2, 4)
    @eval begin
        bitsof(::Type{DNAAlphabet{$n}}) = $n
        bitsof(::Type{RNAAlphabet{$n}}) = $n
    end
end
bitsof(::Type{AminoAcidAlphabet}) = 8
bitsof(::Type{CharAlphabet}) = 32
bitsof(::Type{VoidAlphabet}) = 0

Base.eltype(::Type{DNAAlphabet}) = DNA
Base.eltype(::Type{RNAAlphabet}) = RNA
Base.eltype(::Type{DNAAlphabet{n}}) where {n} = DNA
Base.eltype(::Type{RNAAlphabet{n}}) where {n} = RNA
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid
Base.eltype(::Type{CharAlphabet}) = Char
Base.eltype(::Type{VoidAlphabet}) = Nothing

BioSymbols.alphabet(::Type{DNAAlphabet{2}}) = ACGT
BioSymbols.alphabet(::Type{RNAAlphabet{2}}) = ACGU
BioSymbols.alphabet(::Type{DNAAlphabet{4}}) = alphabet(DNA)
BioSymbols.alphabet(::Type{RNAAlphabet{4}}) = alphabet(RNA)
BioSymbols.alphabet(::Type{AminoAcidAlphabet}) = alphabet(AminoAcid)
# TODO: this alphabet includes invalid Unicode scalar values
BioSymbols.alphabet(::Type{CharAlphabet}) = typemin(Char):typemax(Char)
BioSymbols.alphabet(::Type{VoidAlphabet}) = nothing

# Promotion of Alphabets
# ----------------------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{A}, ::Type{B}) where {A<:$alph,B<:$alph}
        return $alph{max(bitsof(A),bitsof(B))}
    end
end

# Encoders & Decoders
# -------------------

"""
Encode biological characters to binary representation.
"""
function encode end

struct EncodeError{A<:Alphabet,T} <: Exception
    val::T
end

EncodeError(::Type{A}, val::T) where {A,T} = EncodeError{A,T}(val)

function Base.showerror(io::IO, err::EncodeError{A}) where {A}
    print(io, "cannot encode ", err.val, " in ", A)
end

"""
Decode biological characters from binary representation.
"""
function decode end

struct DecodeError{A<:Alphabet,T} <: Exception
    val::T
end

DecodeError(::Type{A}, val::T) where {A,T} = DecodeError{A,T}(val)

function Base.showerror(io::IO, err::DecodeError{A}) where {A}
    print(io, "cannot decode ", err.val, " in ", A)
end


# DNA and RNA alphabets
# ---------------------

for A in (DNAAlphabet, RNAAlphabet)
    T = eltype(A)
    @eval begin
        # 2-bit encoding
        @inline function encode(::Type{$(A){2}}, nt::$(T))
            if count_ones(nt) != 1 || !isvalid(nt)
                throw(EncodeError($(A){2}, nt))
            end
            return convert(UInt8, trailing_zeros(nt))
        end
        @inline function decode(::Type{$(A){2}}, x::UInt8)
            if x > 0x03
                throw(DecodeError($(A){2}, x))
            end
            return reinterpret($(T), 0x01 << x)
        end
        @inline decode(::Type{$(A){2}}, x::Unsigned) = decode($(A){2}, UInt8(x))

        # 4-bit encoding
        @inline function encode(::Type{$(A){4}}, nt::$(T))
            if !isvalid(nt)
                throw(EncodeError($(A){4}, nt))
            end
            return reinterpret(UInt8, nt)
        end
        @inline function decode(::Type{$(A){4}}, x::UInt8)
            if !isvalid($(T), x)
                throw(DecodeError($(A){4}, x))
            end
            return reinterpret($(T), x)
        end
        @inline decode(::Type{$(A){4}}, x::Unsigned) = decode($(A){4}, UInt8(x))
    end
end


# AminoAcidAlphabet
# -----------------

@inline function encode(::Type{AminoAcidAlphabet}, aa::AminoAcid)
    if aa > AA_Gap
        throw(EncodeError(AminoAcidAlphabet, aa))
    end
    return reinterpret(UInt8, aa)
end

@inline function decode(::Type{AminoAcidAlphabet}, x::UInt8)
    if x > 0x1b
        throw(DecodeError(AminoAcidAlphabet, x))
    end
    return reinterpret(AminoAcid, x)
end

@inline function decode(::Type{AminoAcidAlphabet}, x::Unsigned)
    return decode(AminoAcidAlphabet, UInt8(x))
end


# CharAlphabet
# ------------

@inline function encode(::Type{CharAlphabet}, char::Char)
    if char > '\U10ffff'
        throw(EncodeError(CharAlphabet, char))
    end
    return reinterpret(UInt32, char)
end

@inline function decode(::Type{CharAlphabet}, x::UInt32)
    c = reinterpret(Char, x)
    if !isvalid(c)
        throw(DecodeError(CharAlphabet, x))
    end
    return c
end

@inline function decode(::Type{CharAlphabet}, x::Unsigned)
    return decode(CharAlphabet, UInt32(x))
end
