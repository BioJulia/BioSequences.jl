# Alphabet
# ========
#
# Alphabet of biological symbols.
#
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
# Alphabets of biological characters.

An `Alphabet` represents a domain of biological characters.

For example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides
(i.e. A, C, G, and T).

Alphabet types restrict and define the set of biological symbols,
that can be encoded in a given biological sequence type.
They ALSO define *HOW* that encoding is done.

An `Alphabet` type defines the encoding of biological symbols with a pair
of associated `encoder` and `decoder` methods. These paired methods map
between biological symbol values and a binary representation of the symbol.

Any type A <: Alphabet, is expected to implement the `Base.eltype` method
for itself, in addition to a `BioSequences.bitsof` method, and a
`BioSequences.bitsof_t` method. See the docs of `bitsof` and `bitsof_t` for
more detail.

## Required methods and interface
"""
abstract type Alphabet end

"""
Alphabet of nucleic acids.
"""
abstract type NucleicAcidAlphabet{n} <: Alphabet end

"""
DNA nucleotide alphabet.
"""
struct DNAAlphabet{n} <: NucleicAcidAlphabet{n} end

"""
RNA nucleotide alphabet.
"""
struct RNAAlphabet{n} <: NucleicAcidAlphabet{n} end

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

"""
The number of bits required to represent a symbol of the alphabet, in a 
biological sequence.
"""
function bitsof end

bitsof(::Type{A}) where A <: NucleicAcidAlphabet{2} = 2
bitsof(::Type{A}) where A <: NucleicAcidAlphabet{4} = 4
bitsof(::Type{AminoAcidAlphabet}) = 8
bitsof(::Type{CharAlphabet}) = 32
bitsof(::Type{VoidAlphabet}) = 0

"""
The number of bits required to represent a symbol of the alphabet, in a 
biological sequence, as a value type.
"""
function bitsof_t end

bitsof_t(::Type{A}) where A <: NucleicAcidAlphabet{2} = Val{2}
bitsof_t(::Type{A}) where A <: NucleicAcidAlphabet{4} = Val{4}
bitsof_t(::Type{AminoAcidAlphabet}) = Val{8}
bitsof_t(::Type{CharAlphabet}) = Val{32}
bitsof_t(::Type{VoidAlphabet}) = Val{0}

Base.eltype(::Type{A}) where A <: DNAAlphabet = DNA
Base.eltype(::Type{A}) where A <: RNAAlphabet = RNA
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

