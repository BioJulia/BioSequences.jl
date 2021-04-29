```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Sequence Types

BioSequences exports an abstract `BioSequence` type, and several concrete sequence
types which inherit from it.

## The abstract BioSequence

BioSequences provides an abstract type called a `BioSequence{A<:Alphabet}`.
This abstract type, and the methods and traits is supports, allows for
many algorithms in BioSequences to be written as generically as possible,
thus reducing the amount of code to read and understand, whilst maintaining high
performance when such code is compiled for a concrete BioSequence subtype.
Additionally, it allows new types to be implemented that are fully compatible
with the rest of BioSequences, providing that key methods or traits are defined).

This abstract type is parametric over concrete types of `Alphabet`, which
define the range of symbols permitted in the sequence.

Some aliases are also provided for your convenience:

| Type alias      | Type                                 |
| :-------------- | :----------------------------------- |
| `NucleotideSeq` | `BioSequence{<:NucleicAcidAlphabet}` |
| `AminoAcidSeq`  | `BioSequence{AminoAcidAlphabet}`     |

Any concrete sequence type compatible with BioSequences must inherit from
`BioSequence{A}`, where `A` is the alphabet of the concrete sequence type.
It must also have the following methods defined for it:

```@docs
encoded_data
Base.length(::BioSequence)
```

If these requirements are satisfied, the following key traits and methods backing
the BioSequences interface, should be defined already for the sequence type.

```@docs
encoded_data_type
encoded_data_eltype
Alphabet(::BioSequence)
BioSymbols.alphabet(::BioSequence)
BitsPerSymbol
bits_per_symbol
```

As a result, the vast majority of methods described in the rest of this manual
should work out of the box for the concrete sequence type. But they can always
be overloaded if needed.

## Long Sequences

Many genomics scripts and tools benefit from an efficient general purpose
sequence type that allows you to create and edit sequences. In BioSequences,
the `LongSequence` type fills this requirement.

`LongSequence{A<:Alphabet} <: BioSequence{A}` is parameterized by a concrete
`Alphabet` type `A` that defines the domain (or set) of biological symbols
permitted.
For example, `AminoAcidAlphabet` is associated with `AminoAcid` and hence an
object of the `LongSequence{AminoAcidAlphabet}` type represents a sequence of
amino acids.  Symbols from multiple alphabets can't be intermixed in one
sequence type.

The following table summarizes common LongSequence types that have been given
aliases for convenience.

| Type                                | Symbol type | Type alias         |
| :---------------------------------- | :---------- | :----------------- |
| `LongSequence{DNAAlphabet{4}}`      | `DNA`       | `LongDNASeq`       |
| `LongSequence{RNAAlphabet{4}}`      | `RNA`       | `LongRNASeq`       |
| `LongSequence{AminoAcidAlphabet}`   | `AminoAcid` | `LongAminoAcidSeq` |

The `LongDNASeq` and `LongRNASeq` aliases use a DNAAlphabet{4}, which means the
sequence may store ambiguous nucleotides.
If you are sure that nucleotide sequences store unambiguous nucleotides
only, you can reduce the memory required by sequences by using a slightly
different parameter:
`DNAAlphabet{2}` is an alphabet that uses two bits per base and limits to only
unambiguous nucleotide symbols (ACGT in DNA and ACGU in RNA).
Replacing `LongSequence{DNAAlphabet{4}}` in your code with
`LongSequence{DNAAlphabet{2}}` is all you need to do in order to benefit.
Some computations that use bitwise operations will also be dramatically faster.

## Sequence views

Similar to how Base Julia offers views of array objects, BioSequences offers view of
`LongSequence`s - the `LongSubSeq{A<:Alphabet}`.

Conceptually, a `LongSubSeq{A}` is similar to a `LongSequence{A}`, but instead of storing
their own data, they refer to the data of a `LongSequence`. Modiying the `LongSequence`
will be reflected in the view, and vice versa. If the underlying `LongSequence`
is truncated, the behaviour of a view is undefined. For the same reason,
some operations are not supported for views, such as resizing.

The purpose of `LongSubSeq` is that, since they only contain a pointer to the
underlying array, an offset and a length, they are much lighter than `LongSequences`,
and will be stack allocated on Julia 1.5 and newer. Thus, the user may construct
millions of views without major performace implications.

# Alphabet types

```@docs
BioSequences.Alphabet
```

Alphabets control how biological symbols are encoded and decoded.
They also confer many of the automatic traits and methods that any subtype
of `T<:BioSequence{A<:Alphabet}` will get.
