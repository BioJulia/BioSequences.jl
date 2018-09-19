# Alphabet
# ========
#
# Alphabets of biological symbols.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
# Alphabets of biological symbols.

`Alphabet` is perhaps the most important type trait for biological sequences in
BioSequences.jl.

An `Alphabet` represents a domain of biological symbols.

For example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides
(i.e. A, C, G, and T).

Alphabet types restrict and define the set of biological symbols,
that can be encoded in a given biological sequence type.
They ALSO define *HOW* that encoding is done.

An `Alphabet` type defines the encoding of biological symbols with a pair
of associated `encoder` and `decoder` methods. These paired methods map
between biological symbol values and a binary representation of the symbol.

Any type A <: Alphabet, is expected to implement the `Base.eltype` method
for itself.
It is also expected to implement the `BitsPerSymbol` method.
"""
abstract type Alphabet end
Base.eltype(::A) where A <: Alphabet = eltype(A)
Base.firstindex(x::Alphabet) = 1
Base.lastindex(x::Alphabet) = length(x)
symbols(x::Alphabet) = tuple(collect(x)...)


"""
The number of bits required to represent a packed symbol in a vector of bits.
"""
struct BitsPerSymbol{N} end
bits_per_symbol(::BitsPerSymbol{N}) where N = N
bits_per_symbol(::A) where A <: Alphabet = bits_per_symbol(BitsPerSymbol(A()))


# Nucleic acid alphabets
# ====================== 

"""
Alphabet of nucleic acids.
"""
abstract type NucleicAcidAlphabet{n} <: Alphabet end

BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{2} = BitsPerSymbol{2}()
BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{4} = BitsPerSymbol{4}()

@inline function Base.iterate(x::A, state = 0x00) where {A <: NucleicAcidAlphabet{4}}
    state > 0x0F ? nothing : (reinterpret(eltype(x), state), state + 0x01)
end

@inline function Base.iterate(x::A, state = 0x01) where {A <: NucleicAcidAlphabet{2}}
    state > 0x08 ? nothing : (reinterpret(eltype(x), state), state << 0x01)
end

Base.length(x::A) where {A <: NucleicAcidAlphabet{2}} = 4
Base.length(x::A) where {A <: NucleicAcidAlphabet{4}} = 16

@inline function Base.getindex(x::NucleicAcidAlphabet{2}, i::Int)
    return reinterpret(eltype(x), 0x01 << (i - 1))
end
@inline function Base.getindex(x::NucleicAcidAlphabet{4}, i::Int)
    return reinterpret(eltype(x), i - 1)
end


"""
DNA nucleotide alphabet.
"""
struct DNAAlphabet{n} <: NucleicAcidAlphabet{n} end
Base.eltype(::Type{A}) where A <: DNAAlphabet = DNA


"""
RNA nucleotide alphabet.
"""
struct RNAAlphabet{n} <: NucleicAcidAlphabet{n} end
Base.eltype(::Type{A}) where A <: RNAAlphabet = RNA


# Promotion of Alphabets
# ----------------------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{A}, ::Type{B}) where {A<:$alph,B<:$alph}
        # TODO: Resolve this use of bits_per_symbol.
        return $alph{max(bits_per_symbol(A()),bits_per_symbol(B()))}
    end
end


# Amino acid alphabet
# ===================

"""
Amino acid alphabet.
"""
struct AminoAcidAlphabet <: Alphabet end
BitsPerSymbol(::AminoAcidAlphabet) = BitsPerSymbol{8}()
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid
Base.length(x::AminoAcidAlphabet) = 28

@inline function Base.iterate(x::AminoAcidAlphabet, state = 0x00)
    state > 0x1b ? nothing : (reinterpret(AminoAcid, state), state + 0x01)
end

@inline function Base.getindex(x::AminoAcidAlphabet, i)
    return reinterpret(AminoAcid, i - 1)
end


# Generic character alphabet
# ==========================

"""
General character alphabet.
"""
struct CharAlphabet <: Alphabet end
BitsPerSymbol(::CharAlphabet) = BitsPerSymbol{32}()
Base.eltype(::Type{CharAlphabet}) = Char
Base.length(::CharAlphabet) = 1114112

@inline function Base.iterate(::CharAlphabet, state = '\0')
    state > '\U10ffff' ? nothing : (state, state + UInt32(1))
end

@inline function Base.getindex(::CharAlphabet, i)
    return reinterpret(Char, i - 1)
end


# The void alphabet
# =================

"""
Void alphabet (internal use only).
"""
struct VoidAlphabet <: Alphabet end
BitsPerSymbol(::VoidAlphabet) = BitsPerSymbol{0}()
Base.eltype(::Type{VoidAlphabet}) = Nothing
symbols(::VoidAlphabet) = nothing


# Encoders & Decoders
# ===================

"""
Encode biological symbols to binary representation.
"""
function encode end

struct EncodeError{A<:Alphabet,T} <: Exception
    val::T
end

EncodeError(::A, val::T) where {A,T} = EncodeError{A,T}(val)

function Base.showerror(io::IO, err::EncodeError{A}) where {A}
    print(io, "cannot encode ", err.val, " in ", A)
end

"""
Decode biological symbols from binary representation.
"""
function decode end

struct DecodeError{A<:Alphabet,T} <: Exception
    val::T
end

DecodeError(::A, val::T) where {A,T} = DecodeError{A,T}(val)

function Base.showerror(io::IO, err::DecodeError{A}) where {A}
    print(io, "cannot decode ", err.val, " in ", A)
end


# DNA and RNA alphabets
# ---------------------

# A nucleotide with bitvalue B has kmer-bitvalue kmerbits[B+1].
# ambiguous nucleotides have no kmervalue, here set to 0xff
const twobitnucs = (0xff, 0x00, 0x01, 0xff,
                    0x02, 0xff, 0xff, 0xff,
                    0x03, 0xff, 0xff, 0xff,
                    0xff, 0xff, 0xff, 0xff)

for A in (DNAAlphabet, RNAAlphabet)

    T = eltype(A)

    @eval begin

	# 2-bit encoding
        @inline function encode(::$(A){2}, nt::$(T))
            if count_ones(nt) != 1 || !isvalid(nt)
                throw(EncodeError($(A){2}(), nt))
            end
            return convert(UInt8, @inbounds twobitnucs[reinterpret(UInt8, nt) + 0x01])
        end

        @inline function decode(::$(A){2}, x::UInt8)
            if x > 0x03
                throw(DecodeError($(A){2}(), x))
            end
            return reinterpret($(T), 0x01 << x)
        end

	    @inline decode(::$(A){2}, x::Unsigned) = decode($(A){2}(), UInt8(x))

        # 4-bit encoding
        @inline function encode(::$(A){4}, nt::$(T))
            if !isvalid(nt)
                throw(EncodeError($(A){4}(), nt))
            end
            return reinterpret(UInt8, nt)
        end

        @inline function decode(::$(A){4}, x::UInt8)
            if !isvalid($(T), x)
                throw(DecodeError($(A){4}(), x))
            end
            return reinterpret($(T), x)
        end

        @inline decode(::$(A){4}, x::Unsigned) = decode($(A){4}(), UInt8(x))
    end
end


# AminoAcidAlphabet
# -----------------

@inline function encode(::AminoAcidAlphabet, aa::AminoAcid)
    if aa > AA_Gap
        throw(EncodeError(AminoAcidAlphabet(), aa))
    end
    return reinterpret(UInt8, aa)
end

@inline function decode(::AminoAcidAlphabet, x::UInt8)
    if x > 0x1b
        throw(DecodeError(AminoAcidAlphabet(), x))
    end
    return reinterpret(AminoAcid, x)
end

@inline function decode(::AminoAcidAlphabet, x::Unsigned)
    return decode(AminoAcidAlphabet(), UInt8(x))
end


# CharAlphabet
# ------------

@inline function encode(::CharAlphabet, char::Char)
    if char > '\U10ffff'
        throw(EncodeError(CharAlphabet(), char))
    end
    return reinterpret(UInt32, char)
end

@inline function decode(::CharAlphabet, x::UInt32)
    c = reinterpret(Char, x)
    if !isvalid(c)
        throw(DecodeError(CharAlphabet(), x))
    end
    return c
end

@inline function decode(::CharAlphabet, x::Unsigned)
    return decode(CharAlphabet(), UInt32(x))
end
