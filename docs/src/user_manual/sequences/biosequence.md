# Biological sequences

The `BioSequences` module provides representations and tools for manipulating
nucleotide and amino acid sequences.

## Introduction to the sequence data-types


### [The abstract `BioSequence` type](@id biosequence)

BioSequences.jl provides an abstract type called a `BioSequence`.
This abstract type, and the methods and traits is supports, allows for
many algorithms in BioSequences.jl to be created as generic as possible,
thus reducing the amount of code to read and understand, whilst maintaining high
performance when such code is compiled for a concrete BioSequence subtype.
Additionally, it allows new types to be implemented that are fully compatible
with the rest of BioSequences.jl (providing that key methods or traits are defined).

If you are a developer or want to implement a subtype of `BioSequence` you
can find more information on the methods and traits supported by `BioSequence`,
and the requirements of any compliant concrete subtype [here](@ref biosequence_dev).


### Concrete biological sequence types

`BioSequences.jl` currently provides three different concrete sequence types:
`GeneralSequence`, `Kmer` and `ReferenceSequence`.

Different sequence types have different features. In most situations,
`GeneralSequence` type will do and is used as the default representation.
But sometimes other types are much more preferable in terms of memory
efficiency and computation performance.
Here is the summary table of these three types:

| Type                       | Description                                | Element type          | Mutability  | Allocation       |
| :----                      | :-----------                               | :------------         | :---------- | :----------      |
| `BioSequence{A<:Alphabet}` | general-purpose biological sequences       | DNA, RNA, Amino acids | mutable     | heap             |
| `Kmer{T<:NucleicAcid,k}`   | specialized for short nucleotide sequences | DNA, RNA              | immutable   | stack / register |
| `ReferenceSequence`        | specialized for long reference genomes     | DNA                   | immutable   | heap             |

Details specific to these different representations are explained in the following
sections:

* `GeneralSequence`: [General-purpose sequences](@ref)
* `Kmer`: [Nucleic acid k-mers](@ref)
* `ReferenceSequence`: [Reference sequences](@ref)
