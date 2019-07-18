```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# General-purpose sequences

`GeneralSequence{A}` is a generic sequence type parameterized by an
[alphabet](@ref alphabet) type `A` that defines the domain (or set) of
biological symbols, and each alphabet has an associated symbol type.
For example, `AminoAcidAlphabet` is associated with `AminoAcid` and hence an
object of the `GeneralSequence{AminoAcidAlphabet}` type represents a sequence of
amino acids.  Symbols from multiple alphabets can't be intermixed in one
sequence type.

The following table summarizes common sequence types that are defined in the
`BioSequences` module:

| Type                                   | Symbol type | Type alias          |
| :------------------------------------- | :---------- | :------------------ |
| `GeneralSequence{DNAAlphabet{4}}`      | `DNA`       | `LongDNASeq`       |
| `GeneralSequence{RNAAlphabet{4}}`      | `RNA`       | `LongRNASeq`       |
| `GeneralSequence{AminoAcidAlphabet}`   | `AminoAcid` | `LongAminoAcidSeq` |
| `GeneralSequence{CharAlphabet}`        | `Char`      | `LongCharSeq`      |

Parameterized definition of the `GeneralSequence{A}` type is for the purpose of
unifying the data structure and operations of any symbol type. In most cases,
users don't have to care about it and can use *type aliases* listed above.
However, the alphabet type fixes the internal memory encoding and plays an
important role when optimizing performance of a program
(see [Using a more compact sequence representation](@ref) section for low-memory
encodings).  It also enables a developer to define their own alphabet only by
defining few numbers of methods.
This is described in [Defining a new alphabet](@ref) section.


## Constructing `GeneralSequence`s






A random sequence can be obtained by the `randdnaseq`, `randrnaseq` and
`randaaseq` functions, which generate `LongDNASeq`, `LongRNASeq` and
`LongAminoAcidSeq`, respectively. Generated sequences are composed of the
standard symbols without ambiguity and gap. For example, `randdnaseq(6)` may
generate `dna"TCATAG"` but never generates `dna"TNANAG"` or `dna"T-ATAG"`.

A translatable `LongRNASeq` can also be converted to an `LongAminoAcidSeq`
using the [`translate`](@ref) function.






### Setindex and modifying DNA sequences





## Site counting

BioSequences extends the `Base.count` method to provide some useful utilities for
counting the number of sites in biological sequences.

### Site types

Different types of site can be counted. Each of these types is a concrete
subtype of the abstract type `Site`:

```@docs
Certain
Gap
Ambiguous
Match
Mismatch
```

### `Base.count` methods

The count method can be used with two sequences and a concrete subtype of
`Site`:

```jldoctest
julia> count(Match, dna"ATCGATCG", dna"AAGGTTCG")
5
```

By providing a `window` and `step` size, counting can be done from within
a sliding window:

```jldoctest
julia> count(Match, dna"ATCGATCG", dna"AAGGTTCG", 3, 1)
6-element Array{IntervalTrees.IntervalValue{Int64,Int64},1}:
 IntervalTrees.IntervalValue{Int64,Int64}
(1,3) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(2,4) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(3,5) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(4,6) => 2
 IntervalTrees.IntervalValue{Int64,Int64}
(5,7) => 2
 IntervalTrees.IntervalValue{Int64,Int64}
(6,8) => 3
```

### The `pairwise_count` function
Counting can also be done on a set of sequences in a pairwise manner with the
`count_pairwise` function:

```jldoctest
julia> count_pairwise(Match, dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA")
4×4 Array{Int64,2}:
 0  6  7  5
 6  0  7  7
 7  7  0  6
 5  7  6  0
```

## Iteration

Sequences also work as iterators over symbols:

```jldoctest
julia> n = 0
0

julia> for nt in dna"ATNGNNT"
           if nt == DNA_N
               n += 1
           end
       end

julia> n
3

```




Recall that `LongDNASeq` is a type alias of `GeneralSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `GeneralSequence{DNAAlphabet{2}}` saves half
of memory footprint compared to `GeneralSequence{DNAAlphabet{4}}`. If you need to
handle reference genomes that are composed of five nucleotides, ACGTN,
consider to use the `ReferenceSequence` type described in the [Reference
sequences](@ref) section.
