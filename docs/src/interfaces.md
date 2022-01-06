```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Custom BioSequences types

If you're a developing your own Bioinformatics package or method, you may find
that the reference implementation of concrete `LongSequence` types provided in
this package are not optimal for your purposes.

This page describes the interfaces for BioSequences' core types for
developers or other packages implementing their own sequence types or extending
BioSequences functionality.

## Implementing custom Alphabets

Recall the required methods that define the [`Alphabet`](@ref) interface. 

To create an example custom alphabet, we need to create a singleton type, that
implements a few methods in order to conform to the interface as described in the
[`Alphabet`](@ref) documentation.

Let's do that for a restricted Amino Acid alphabet. We can test that it conforms
to the interface with the [`BioSequences.has_interface`](@ref) function.

```jldoctest
julia> struct ReducedAAAlphabet <: Alphabet end

julia> Base.eltype(::Type{ReducedAAAlphabet}) = AminoAcid

julia> BioSequences.BitsPerSymbol(::ReducedAAAlphabet) = BioSequences.BitsPerSymbol{4}()

julia> function BioSequences.symbols(::ReducedAAAlphabet)
           (AA_L, AA_C, AA_A, AA_G, AA_S, AA_T, AA_P, AA_F,
            AA_W, AA_E, AA_D, AA_N, AA_Q, AA_K, AA_H, AA_M)
       end

julia> const (ENC_LUT, DEC_LUT) = let
           enc_lut = fill(0xff, length(alphabet(AminoAcid)))
           dec_lut = fill(AA_A, length(symbols(ReducedAAAlphabet())))
           for (i, aa) in enumerate(symbols(ReducedAAAlphabet()))
               enc_lut[reinterpret(UInt8, aa) + 0x01] = i - 1
               dec_lut[i] = aa
           end
           (Tuple(enc_lut), Tuple(dec_lut))
       end
((0x02, 0xff, 0x0b, 0x0a, 0x01, 0x0c, 0x09, 0x03, 0x0e, 0xff, 0x00, 0x0d, 0x0f, 0x07, 0x06, 0x04, 0x05, 0x08, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff), (AA_L, AA_C, AA_A, AA_G, AA_S, AA_T, AA_P, AA_F, AA_W, AA_E, AA_D, AA_N, AA_Q, AA_K, AA_H, AA_M))

julia> function BioSequences.encode(::ReducedAAAlphabet, aa::AminoAcid)
           i = reinterpret(UInt8, aa) + 0x01
           (i ≥ length(ENC_LUT) || @inbounds ENC_LUT[i] === 0xff) && throw(DomainError(aa))
           (@inbounds ENC_LUT[i]) % UInt
       end

julia> function BioSequences.decode(::ReducedAAAlphabet, x::UInt)
           x ≥ length(DEC_LUT) && throw(DomainError(aa))
           @inbounds DEC_LUT[x + UInt(1)]
       end

julia> BioSequences.has_interface(Alphabet, ReducedAAAlphabet())
true

```

## Implementing custom BioSequences

Recall the required methods that define the [`BioSequence`](@ref) interface. 

To create an example custom alphabet, we need to create a singleton type, that
implements a few methods in order to conform to the interface as described in the
[`BioSequence`](@ref) documentation.

Let's do that for a custom sequence type that is optimised to represent a small
sequence: A Codon. We can test that it conforms to the interface with the
[`BioSequences.has_interface`](@ref) function.

```jldoctest
julia> struct Codon <: BioSequence{RNAAlphabet{2}}
           x::UInt8
       end

julia> function Codon(iterable)
           length(iterable) == 3 || error("Must have length 3")
           x = zero(UInt)
           for (i, nt) in enumerate(iterable)
               x |= BioSequences.encode(Alphabet(Codon), convert(RNA, nt)) << (6-2i)
           end
           Codon(x % UInt8)
       end
Codon

julia> Base.length(::Codon) = 3

julia> BioSequences.encoded_data_eltype(::Type{Codon}) = UInt

julia> function BioSequences.extract_encoded_element(x::Codon, i::Int)
           ((x.x >>> (6-2i)) & 3) % UInt
       end

julia> Base.copy(seq::Codon) = Codon(seq.x)
       
julia> BioSequences.has_interface(BioSequence, Codon, [RNA_C, RNA_U, RNA_A], false)
true
```

## Interface checking functions

```@docs
BioSequences.has_interface
```