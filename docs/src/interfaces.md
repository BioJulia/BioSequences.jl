# BioSequences Interfaces

This page describes the interfaces for BioSequences' core types for
developers or other packages implementing their own sequence types or extending
BioSequences functionality.

## Alphabets

`Alphabet` is an abstract type defined by BioSequences.

Subtypes of `Alphabet` are singleton structs that may or may not be parameterized.

Alphabets determine a _finite_ list of biological symbols that can be encoded in
the given sequence type.

`Alphabet`s controls the encoding/decoding between a decoded element type and an
internal data representation type.

An `Alphabet` must never encode (using `encode`) or decode (using `decode`) invalid
data, but must error in cases where invalid data could have been produced.
Other methods for check-free encoding/decoding methods may be added.

### Required implementation

Every subtype `A` of `Alphabet` must implement:

- `Base.eltype(::Type{A})::Type{E}` for some eltype `E`.
- `symbols(::A)::Tuple{Vararg{E}}`. This gives an ordered tuples of all elements of `A`.
- `encode(::A, ::E)::X` encodes an element to the internal data eltype `X`
- `decode(::A, ::X)::E` decodes an `X` to an element `E`.
- Except for `eltype` which must follow Base conventions, all functions
  operating on Alphabet should operate on instances of the alphabet, not the type.

If you want interoperation with existing subtypes of `BioSequence`, the element type E must be UInt, and you must also implement:

- `BitsPerSymbol(::A)::BitsPerSymbol{N}`, where the `N` must be zero or a power
  of two in `[1, 8 * sizeof(UInt)]`.

#### Minimal example

This example implements a restricted AminoAcid alphabet, that does not contain
_all_ possible IUPAC symbols.

```julia
struct ReducedAAAlphabet <: Alphabet end

Base.eltype(::Type{ReducedAAAlphabet}) = AminoAcid

BioSequences.BitsPerSymbol(::ReducedAAAlphabet) = BioSequences.BitsPerSymbol{4}()

function BioSequences.symbols(::ReducedAAAlphabet)
    (AA_L, AA_C, AA_A, AA_G, AA_S, AA_T, AA_P, AA_F,
     AA_W, AA_E, AA_D, AA_N, AA_Q, AA_K, AA_H, AA_M)
end

const (ENC_LUT, DEC_LUT) = let
    enc_lut = fill(0xff, length(alphabet(AminoAcid)))
    dec_lut = fill(AA_A, length(symbols(ReducedAAAlphabet())))
    for (i, aa) in enumerate(symbols(ReducedAAAlphabet()))
        enc_lut[reinterpret(UInt8, aa) + 0x01] = i - 1
        dec_lut[i] = aa
    end
    (Tuple(enc_lut), Tuple(dec_lut))
end

function BioSequences.encode(::ReducedAAAlphabet, aa::AminoAcid)
    i = reinterpret(UInt8, aa) + 0x01
    (i ≥ length(ENC_LUT) || @inbounds ENC_LUT[i] === 0xff) && throw(DomainError(aa))
    (@inbounds ENC_LUT[i]) % UInt
end

function BioSequences.decode(::ReducedAAAlphabet, x::UInt)
    x ≥ length(DEC_LUT) && throw(DomainError(aa))
    @inbounds DEC_LUT[x + UInt(1)]
end
```

## BioSequences

`BioSequence` subtypes are linear container types, with random access and indices
`Base.OneTo(length(x))`.

They contain zero or more coding elements of type `encoded_data_eltype(typeof(x))`.

They may or may not be mutable, therefore generic functions cannot assume mutability.

They are associated with an `Alphabet`, `A` by being a subtype of `BioSequence{A}`.

The concrete type, together with its associated `Alphabet`, and no other
properties, determines how to extract individual elements from the encoded data
and the index.

The `BioSequence` subtype with the index, and optionally with the `Alphabet`,
determines how to extract the internally encoded elements.

The `Alphabet` then controls how to decode it to an element.

### Required implementation

Every subtype `T` of `BioSequence{A}` must implement the following methods,
where x in an instance of the subtype:

- `Base.length(x)::Int`
- `encoded_data_eltype(::Type{T})::Type{E}`
- `extract_encoded_element(x, ::Int)::E`
- `T` must be able to be constructed from any iterable with length defined and
  with a known, compatible element type.

Furthermore, _mutable_ sequences should implement

- `encoded_setindex!(x, ::E, ::Int)`

For compatibility with existing Alphabets, the encoded data eltype must be `UInt`.

#### Minimal example

```julia
struct Codon <: BioSequence{RNAAlphabet{2}}
    x::UInt8
end

function Codon(iterable)
    length(iterable) == 3 || error("Must have length 3")
    x = zero(UInt)
    for (i, nt) in enumerate(iterable)
        x |= BioSequences.encode(Alphabet(Codon), convert(RNA, nt)) << (6-2i)
    end
    Codon(x)
end 

Base.length(::Codon) = 3
BioSequences.encoded_data_eltype(::Type{Codon}) = UInt
function BioSequences.extract_encoded_element(x::Codon, i::Int)
    (x.x >>> ((i-1) & 7)) % UInt
end
```