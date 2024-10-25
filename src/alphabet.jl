###
### Alphabet
###
###
### Alphabets of biological symbols.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    Alphabet

`Alphabet` is the most important type trait for `BioSequence`. An `Alphabet`
represents a set of biological symbols encoded by a sequence, e.g. A, C, G
and T for a DNA Alphabet that requires only 2 bits to represent each symbol.

# Extended help
* Subtypes of Alphabet are singleton structs that may or may not be parameterized.
* Alphabets span over a *finite* set of biological symbols.
* The alphabet controls the encoding from some internal "encoded data" to a BioSymbol 
  of the alphabet's element type, as well as the decoding, the inverse process.
* An `Alphabet`'s `encode` method must not produce invalid data. 

### Required methods
Every subtype `A` of `Alphabet` must implement:
* `Base.eltype(::Type{A})::Type{S}` for some eltype `S`, which must be a `BioSymbol`.
* `symbols(::A)::Tuple{Vararg{S}}`. This gives tuples of all symbols in the set of `A`.
* `encode(::A, ::S)::E` encodes a symbol to an internal data eltype `E`.
* `decode(::A, ::E)::S` decodes an internal data eltype `E` to a symbol `S`.
* Except for `eltype` which must follow Base conventions, all functions operating
  on `Alphabet` should operate on instances of the alphabet, not the type.

If you want interoperation with existing subtypes of `BioSequence`,
the encoded representation `E` must be of type `UInt`, and you must also implement:
* `BitsPerSymbol(::A)::BitsPerSymbol{N}`, where the `N` must be zero
  or a power of two in [1, 2, 4, 8, 16, 32, [64 for 64-bit systems]].

### Optional methods
* `BitsPerSymbol` for compatibility with existing `BioSequence`s
* `AsciiAlphabet` for increased printing/writing efficiency
* `tryencode` for fallible encoding.
"""
abstract type Alphabet end

"""
    function has_interface(::Type{Alphabet}, A::Alphabet)

Returns whether `A` conforms to the `Alphabet` interface.
"""
function has_interface(::Type{Alphabet}, A::Alphabet)
    try
        eltype(typeof(A)) <: BioSymbol || return false
        syms = symbols(A)
        (syms isa Tuple{Vararg{eltype(typeof(A))}} && length(syms) > 0) || return false
        codes = map(i -> encode(A, i), syms)
        codes isa (NTuple{N, T} where {N, T <: Union{UInt8, UInt16, UInt32, UInt64}}) || return false
        recodes = map(i -> decode(A, i), codes)
        syms == recodes || return false
        bps = BitsPerSymbol(A)
        bps isa BitsPerSymbol || return false
        in(BioSequences.bits_per_symbol(A), (0, 1, 2, 4, 8, 16, 32, 64)) || return false
    catch error
        error isa MethodError && return false
        rethrow(error)
    end
    return true
end

"""
The number of bits required to represent a packed symbol encoding in a vector of bits.
"""
bits_per_symbol(A::Alphabet) = bits_per_symbol(BitsPerSymbol(A))
Base.length(A::Alphabet) = length(symbols(A))

"""
    BitsPerSymbol{N}
A trait object specifying the number of bits it takes to encode a biosymbol in an `Alphabet`
Alphabets `A` should implement `BitsPerSymbol(::A)`.
For compatibility with existing BioSequences, the number of bits should be a power of two
between 1 and 32, both inclusive.
See also: [`Alphabet`](@ref)
"""
struct BitsPerSymbol{N} end
bits_per_symbol(::BitsPerSymbol{N}) where N = N

"Compute whether all bitpatterns represent valid symbols for an alphabet"
iscomplete(A::Alphabet) = Val(length(symbols(A)) === 1 << bits_per_symbol(A))

## Encoders & Decoders

"""
    encode(::Alphabet, s::BioSymbol)

Internal function, do not use in user code.
Encode BioSymbol `s` to an internal representation using an [`Alphabet`](@ref).
This decoding is checked to enforce valid data element.
If `s` cannot be encoded to the given alphabet, throw an `EncodeError`

"""
@inline function encode(A::Alphabet, s::BioSymbol)
    y = @inline tryencode(A, s)
    return y === nothing ? throw(EncodeError(A, s)) : y
end

tryencode(A::Alphabet, s::Any) = nothing

"""
    tryencode(::Alphabet, x::S)
Try encoding BioSymbol `S` to the internal representation of [`Alphabet`](@ref),
returning `nothing` if not successful.

See also: `encode`[@ref], `decode`[@ref]
"""
function tryencode end

"""
    EncodeError

Exception thrown when a `BioSymbol` cannot be encoded to a given [`Alphabet`](@ref).

# Examples
```
julia> try
           BioSequences.encode(DNAAlphabet{2}(), DNA_N)
       catch err
           println(err isa BioSequences.EncodeError)
       end
true
```
"""
struct EncodeError{A<:Alphabet,T} <: Exception
    val::T
end

EncodeError(::A, val::T) where {A,T} = EncodeError{A,T}(val)

function Base.showerror(io::IO, err::EncodeError{A}) where {A}
    print(io, "cannot encode ", repr(err.val), " in ", A)
end

"""
    decode(::Alphabet, x::E)

Decode internal representation `E` to a `BioSymbol` using an `Alphabet`.
"""
function decode end

function Base.iterate(a::Alphabet, state = 1)
	state > length(a) && return nothing
	@inbounds sym = symbols(a)[state]
	return sym, state + 1
end

## Nucleic acid alphabets

"""
Alphabet of nucleic acids. Parameterized by the number of bits per symbol, by
default only `2` or `4`-bit variants exists.
"""
abstract type NucleicAcidAlphabet{N} <: Alphabet end

"""
DNA nucleotide alphabet.

`DNAAlphabet` has a parameter `N` which is a number that determines the
`BitsPerSymbol` trait. Currently supported values of `N` are 2 and 4.
"""
struct DNAAlphabet{N} <: NucleicAcidAlphabet{N} end
Base.eltype(::Type{<:DNAAlphabet}) = DNA

"""
RNA nucleotide alphabet.

`RNAAlphabet` has a parameter `N` which is a number that determines the
`BitsPerSymbol` trait. Currently supported values of `N` are 2 and 4.
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


## Encoding and decoding DNA and RNA alphabet symbols
for A in (DNAAlphabet, RNAAlphabet)
    T = eltype(A)
    @eval begin

	    # 2-bit encoding
        function tryencode(::$(A){2}, nt::$(T))
            u = reinterpret(UInt8, nt)
            isone(count_ones(u)) || return nothing
            trailing_zeros(u) % UInt
        end

        @inline function decode(::$(A){2}, x::UInt)
            return reinterpret($(T), 0x01 << (x & 0x03))
        end

	    @inline decode(::$(A){2}, x::Unsigned) = decode($(A){2}(), UInt(x))

        # 4-bit encoding
        function tryencode(::$(A){4}, nt::$(T))
            return convert(UInt, reinterpret(UInt8, nt))
        end

        @inline function decode(::$(A){4}, x::UInt)
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

function tryencode(::AminoAcidAlphabet, aa::AminoAcid)
    reinterpret(UInt8, aa) > reinterpret(UInt8, AA_Gap) && return nothing
    return convert(UInt, reinterpret(UInt8, aa))
end

@inline function decode(::AminoAcidAlphabet, x::UInt)
    return reinterpret(AminoAcid, x % UInt8)
end

@inline function decode(::AminoAcidAlphabet, x::Unsigned)
    return decode(AminoAcidAlphabet(), UInt(x))
end

# AsciiAlphabet trait - add to user defined type to use speedups.
# Must define methods codetype, BioSymbols.stringbyte, ascii_encode
"Abstract trait for ASCII/Unicode dispatch. See `AsciiAlphabet`"
abstract type AlphabetCode end

"""
    AsciiAlphabet

Trait for alphabet using ASCII characters as String representation.
Define `codetype(A) = AsciiAlphabet()` for a user-defined `Alphabet` A to gain speed.
Methods needed: `BioSymbols.stringbyte(::eltype(A))` and `ascii_encode(A, ::UInt8)`.
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

"""
	ascii_encode(::Alphabet, b::UInt8)::UInt8

Encode the ASCII character represented by `b` to the internal alphabet encoding.
For example, the input byte `UInt8('C')` is encoded to `0x01` and `0x02` for
2- and 4-bit DNA alphabets, reprectively.
This method is only needed if the `Alphabet` is an `AsciiAlphabet`.

See also: [`BioSequences.AsciiAlphabet`](@ref)
"""
function ascii_encode end

for (anum, atype) in enumerate((DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4},
    RNAAlphabet{2}, AminoAcidAlphabet))
    tablename = Symbol("BYTE_TO_ALPHABET_CHAR" * string(anum))
    @eval begin
        alph = $(atype)()
        syms = symbols(alph)
        const $(tablename) = let
            bytes = fill(0x80, 256)
            for symbol in syms
                bytes[UInt8(Char(symbol)) + 1] = encode(alph, symbol)
                bytes[UInt8(lowercase(Char(symbol))) + 1] = encode(alph, symbol)
            end
            Tuple(bytes)
        end
        ascii_encode(::$(atype), x::UInt8) = @inbounds $(tablename)[x + 1]
    end
end

const GUESS_ALPHABET_LUT = let
    v = zeros(UInt8, 64)
    for (offset, A) in [
        (0, DNAAlphabet{2}()),
        (0, RNAAlphabet{2}()),
        (1, DNAAlphabet{4}()),
        (1, RNAAlphabet{4}()),
        (2, DNAAlphabet{4}()),
        (3, AminoAcidAlphabet())
    ]
        for symbol in A
            for byte in (UInt8(uppercase(Char(symbol))), UInt8(lowercase(Char(symbol))))
                v[div(byte, 2) + 1] |= 0x01 << (4*isodd(byte) + offset)
            end
        end
    end
    Tuple(v)
end

"""
    guess_alphabet(s::Union{AbstractString, AbstractVector{UInt8}}) -> Union{Integer, Alphabet}

Pick an `Alphabet` that can encode input `s`.  If no `Alphabet` can, return the index of the first
byte of the input which is not encodable in any alphabet.
This function only knows about the alphabets listed below. If multiple alphabets are possible,
pick the first from the order below (i.e. `DNAAlphabet{2}()` if possible, otherwise `RNAAlphabet{2}()` etc).
1. `DNAAlphabet{2}()`
2. `RNAAlphabet{2}()`
3. `DNAAlphabet{4}()`
4. `RNAAlphabet{4}()`
5. `AminoAcidAlphabet()`

!!! warning
    The functions `bioseq` and `guess_alphabet` are intended for use in interactive
    sessions, and are not suitable for use in packages or non-ephemeral work.
    They are type unstable, and their heuristics **are subject to change** in minor versions.

# Examples
```jldoctest
julia> guess_alphabet("AGGCA")
DNAAlphabet{2}()

julia> guess_alphabet("WKLQSTV")
AminoAcidAlphabet()

julia> guess_alphabet("QAWT+!")
5

julia> guess_alphabet("UAGCSKMU")
RNAAlphabet{4}()
```
"""
function guess_alphabet(v::AbstractVector{UInt8})
    possibilities = 0x0f
    for (index, byte) in pairs(v)
        lut_byte = @inbounds GUESS_ALPHABET_LUT[div(byte & 0x7f, 2) + 0x01]
        lut_value = (lut_byte >>> (4 * isodd(byte))) & 0x0f
        possibilities &= (lut_value * (byte < 0x80))
        iszero(possibilities) && return index
    end
    dna = !iszero(possibilities & 0b0100)
    @assert !iszero(possibilities) # We checked that in the loop above
    if !iszero(possibilities & 0b0001)
        dna ? DNAAlphabet{2}() : RNAAlphabet{2}()
    elseif !iszero(possibilities & 0b0010)
        dna ? DNAAlphabet{4}() : RNAAlphabet{4}()
    else
        AminoAcidAlphabet()
    end
end
guess_alphabet(s::AbstractString) = guess_alphabet(codeunits(s))

const KNOWN_ALPHABETS = Union{DNAAlphabet, RNAAlphabet, AminoAcidAlphabet}
