###
### Alphabet
###
###
### Alphabets of biological symbols.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
# Alphabets of biological symbols.

`Alphabet` is the most important type trait for `BioSequence` An `Alphabet`
represents a set of biological symbols encoded by a sequence, e.g. A, C, G
and T for a 2-bit DNA Alphabet.

* Subtypes of Alphabet are singleton structs that may or may not be parameterized.
* Alphabets span over a *finite* list of biological symbols
* The alphabet controls the encoding/decoding between a decoded element type and
    an internal data representation type.
* An `Alphabet` must never encode (using `encode`) or decode (using `decode`)
invalid data.  Other methods for check-free encoding/decoding methods may be added.

Every subtype `A` of `Alphabet` must implement:
* `Base.eltype(::Type{A})::Type{E}` for some eltype `E`, which must be a BioSymbol
* `symbols(::A)::Tuple{Vararg{E}}`. This gives tuplea of all elements of `A`.
* `encode(::A, ::E)::X` encodes an element to the internal data eltype `X`
* `decode(::A, ::X)::E` decodes an `X` to an element `E`.
* Except for `eltype` which must follow Base conventions, all functions operating
on `Alphabet` should operate on instances of the alphabet, not the type.


If you want interoperation with existing subtypes of `BioSequence`,
the element type `E` must be `UInt`, and you must also implement:
* `BitsPerSymbol(::A)::BitsPerSymbol{N}`, where the `N` must be zero
or a power of two in [1, 2, 4, 8, 16, 32, [64 for 64-bit systems]].
"""
abstract type Alphabet end

"""
The number of bits required to represent a packed symbol in a vector of bits.
"""
bits_per_symbol(A::Alphabet) = bits_per_symbol(BitsPerSymbol(A))
Base.length(A::Alphabet) = length(symbols(A))

## Bits per symbol

struct BitsPerSymbol{N} end
bits_per_symbol(::BitsPerSymbol{N}) where N = N

## Encoders & Decoders

"""
    encode(::Alphabet, x::T)

Encode biosymbol `T` to a data element using an `Alphabet`.
This decoding is checked to enforce valid data element.
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
    decode(::Alphabet, x::E)

Decode data element `E` to a `BioSymbol` using an `Alphabet`.
This decoding is checked to enforce valid biosymbols.
"""
function decode end

struct DecodeError{A<:Alphabet,T} <: Exception
    val::T
end

DecodeError(::A, val::T) where {A,T} = DecodeError{A,T}(val)

function Base.showerror(io::IO, err::DecodeError{A}) where {A}
    print(io, "cannot decode ", err.val, " in ", A)
end

## Nucleic acid alphabets

"""
Alphabet of nucleic acids.
"""
abstract type NucleicAcidAlphabet{N} <: Alphabet end

"""
DNA nucleotide alphabet.
"""
struct DNAAlphabet{N} <: NucleicAcidAlphabet{N} end
Base.eltype(::Type{<:DNAAlphabet}) = DNA

"""
RNA nucleotide alphabet.
"""
struct RNAAlphabet{N} <: NucleicAcidAlphabet{N} end
Base.eltype(::Type{<:RNAAlphabet}) = RNA

symbols(::DNAAlphabet{2}) = (DNA_A, DNA_C, DNA_G, DNA_T)
symbols(::RNAAlphabet{2}) = (RNA_A, RNA_C, RNA_G, RNA_U)

function symbols(::DNAAlphabet{4})
    (DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V,
    DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N)
end

function symbols(::RNAAlphabet{4})
    (RNA_Gap, RNA_A, RNA_C, RNA_M, RNA_G, RNA_R, RNA_S, RNA_V,
    RNA_U, RNA_W, RNA_Y, RNA_H, RNA_K, RNA_D, RNA_B, RNA_N)
end

BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{2} = BitsPerSymbol{2}()
BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{4} = BitsPerSymbol{4}()


## Encoding and decoding DNA and RNA alphabets

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
            return convert(UInt, @inbounds twobitnucs[reinterpret(UInt8, nt) + 0x01])
        end

        @inline function decode(::$(A){2}, x::UInt)
            if x > UInt(3)
                throw(DecodeError($(A){2}(), x))
            end
            return reinterpret($(T), 0x01 << (x & 0x03))
        end

	    @inline decode(::$(A){2}, x::Unsigned) = decode($(A){2}(), UInt(x))

        # 4-bit encoding
        @inline function encode(::$(A){4}, nt::$(T))
            if !isvalid(nt)
                throw(EncodeError($(A){4}(), nt))
            end
            return convert(UInt, reinterpret(UInt8, nt))
        end

        @inline function decode(::$(A){4}, x::UInt)
            if !isvalid($(T), x)
                throw(DecodeError($(A){4}(), x))
            end
            return reinterpret($(T), x % UInt8)
        end

        @inline decode(::$(A){4}, x::Unsigned) = decode($(A){4}(), UInt(x))
    end
end

### Promotion of nucletid acid alphabets

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{A}, ::Type{B}) where {A<:$alph,B<:$alph}
        # TODO: Resolve this use of bits_per_symbol.
        return $alph{max(bits_per_symbol(A()),bits_per_symbol(B()))}
    end
end

## Amino acid alphabet

"""
Amino acid alphabet.
"""
struct AminoAcidAlphabet <: Alphabet end
BitsPerSymbol(::AminoAcidAlphabet) = BitsPerSymbol{8}()
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid

function symbols(::AminoAcidAlphabet)
    (AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H,
    AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W,
    AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap)
end

@inline function encode(::AminoAcidAlphabet, aa::AminoAcid)
    if reinterpret(UInt8, aa) > reinterpret(UInt8, AA_Gap)
        throw(EncodeError(AminoAcidAlphabet(), aa))
    end
    return convert(UInt, reinterpret(UInt8, aa))
end

@inline function decode(::AminoAcidAlphabet, x::UInt)
    if x > 0x1b
        throw(DecodeError(AminoAcidAlphabet(), x))
    end
    return reinterpret(AminoAcid, x % UInt8)
end

@inline function decode(::AminoAcidAlphabet, x::Unsigned)
    return decode(AminoAcidAlphabet(), UInt(x))
end

# AsciiAlphabet trait - add to user defined type to use speedups.
# Must define methods codetype, stringbyte,
"Abstract trait for ASCII/Unicode dispatch. See `AsciiAlphabet`"
abstract type AlphabetCode end

"""
    AsciiAlphabet

Trait for alphabet using ASCII characters as String representation.
Define `codetype(A) = AsciiAlphabet()` for a user-defined `Alphabet` A to gain speed.
Methods needed: `stringbyte(::eltype(A))` and `stringbyte(A, ::UInt8)`.
"""
struct AsciiAlphabet <: AlphabetCode end

"Trait for alphabet using Unicode. See `AsciiAlphabet`"
struct UnicodeAlphabet <: AlphabetCode end

function codetype(::A) where {A <: Union{DNAAlphabet{2}, DNAAlphabet{4},
                                         RNAAlphabet{2}, RNAAlphabet{4},
                                         AminoAcidAlphabet}}
    return AsciiAlphabet()
end
codetype(::Alphabet) = UnicodeAlphabet()

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
