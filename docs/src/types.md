```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Abstract Types

BioSequences exports an abstract `BioSequence` type and several concrete sequence
types that inherit from it.

## The abstract BioSequence

BioSequences provides an abstract type called a `BioSequence{A<:Alphabet}`.
This abstract type, and the methods and traits it supports, allow for
many algorithms in BioSequences to be written as generically as possible,
thus reducing the amount of code to read and understand, whilst maintaining high
performance when such code is compiled for a concrete BioSequence subtype.
Additionally, it allows new types to be implemented that are fully compatible
with the rest of BioSequences, provided that key methods or traits are defined.

```@docs
BioSequence
```

Some aliases for `BioSequence` are also provided for your convenience:

```@docs
NucSeq
AASeq
```

Let's have a closer look at some of the methods that a subtype of `BioSequence`
must implement. Check out the Julia Base library docs for `length`, `copy`, and `resize!`.

```@docs
encoded_data_eltype
extract_encoded_element
encoded_setindex!
```

A correctly defined subtype of `BioSequence` that satisfies the interface should
have the vast majority of methods described in the rest of this manual work out
of the box for that type. But they can always be overloaded if needed. Indeed,
the `LongSequence` type overloads some of the generic `BioSequence` methods, for
example transformation and counting operations where efficiency gains can be
made due to the specific internal representation of the type.


## The abstract Alphabet

Alphabets control how biological symbols are encoded and decoded.
They also confer many of the automatic traits and methods that any subtype
of `T<:BioSequence{A<:Alphabet}` will get.

```@docs
BioSequences.Alphabet
BioSequences.AsciiAlphabet
```

# Concrete types

## Implemented alphabets

```@docs
DNAAlphabet
RNAAlphabet
AminoAcidAlphabet
```

## Long Sequences

```@docs
LongSequence
```

## Sequence views

Similar to how Julia's `Base` offers views of array objects, BioSequences offers views of
`LongSequence`s—the `LongSubSeq{A<:Alphabet}`.

Conceptually, a `LongSubSeq{A}` is similar to a `LongSequence{A}`, but instead of storing
its own data, it refers to the data of a `LongSequence`. Modifying the `LongSequence`
will be reflected in the view, and vice versa. If the underlying `LongSequence`
is truncated, the behaviour of a view is undefined. For the same reason,
some operations are not supported for views, such as resizing.

The purpose of `LongSubSeq` is that, since it only contains a pointer to the
underlying array, an offset, and a length, it is much lighter than `LongSequence`s,
and will be stack-allocated on Julia 1.5 and newer. Thus, the user may construct
millions of views without major performance implications.
