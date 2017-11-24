# Biological sequences

The `BioSequences` module provides representations and tools for manipulating
nucleotide and amino acid sequences.

## Introduction to the sequence data-types


### The abstract `Sequence` type

BioSequences.jl provides an abstract type called a `Sequence`.
This abstract type allows for some in BioSequences.jl to be created as generic
as possible, thus reducing the amount of code to read and understand, and also
allowing new types to be concieved that are fully compatible with the rest of
BioSequences.jl (providing that key methods or traits are defined).

docs@```
Sequence
```

In addition to being compatible with many `Base` methods.
Any `Sequence` supports several generalised methods.
These methods are sometimes overloaded (redefined) for the concrete
sequence types, if some internal aspect of their representation allows
for a much more performant version of the method.
Some of these generic methods also make use of some of the traits `Sequences`
support in order to specialise their behaviour.

docs@```
gc_content
count_gc
ungap
ungap!
isrepetetive
ispalindromic
hasambiguity
```

### Concrete biological sequence types

Sequences in BioSequences.jl are more strictly typed than in many other
libraries;
Elements in a sequence are typed as biological symbol instead of character or
byte. They are special-purpose types rather than simply strings and hence
offer additional functionality that naive string types don't have.

Though this strictness sacrifices some convenience, it also means you can
always rely on a DNA sequence type to store DNA and nothing but DNA, without
having to check, or deal with lowercase versus uppercase and so on.

Strict separation of sequence types also means we are free to choose the most
efficient representation. DNA and RNA sequences are encoded using either four
bits per base (which is the default), or two bits per base. This makes them
memory efficient and allows us to speed up many common operations and
transformations, like nucleotide composition, reverse complement, and *k*-mer
enumeration.

`BioSequences.jl` provides three different concrete sequence types:
`BioSequence`, `Kmer` and `ReferenceSequence`.

Different sequence types have different features. In most situations,
`BioSequence` type will do and is used as the default representation.
But sometimes other types are much more preferable in terms of memory
efficiency and computation performance.
Here is the summary table of these three types:

| Type                       | Description                                | Element type          | Mutability  | Allocation       |
| :----                      | :-----------                               | :------------         | :---------- | :----------      |
| `BioSequence{A<:Alphabet}` | general-purpose biological sequences       | DNA, RNA, Amino acids | mutable     | heap             |
| `Kmer{T<:NucleicAcid,k}`   | specialized for short nucleotide sequences | DNA, RNA              | immutable   | stack / register |
| `ReferenceSequence`        | specialized for long reference genomes     | DNA                   | immutable   | heap             |

Details of these different representations are explained in the following
sections:

* `BioSequence`: [General-purpose sequences](@ref)
* `Kmer`: [Nucleic acid k-mers](@ref)
* `ReferenceSequence`: [Reference sequences](@ref)
