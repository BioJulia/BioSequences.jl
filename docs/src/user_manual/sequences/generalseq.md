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

```
julia> count(Match, dna"ATCGATCG", dna"AAGGTTCG")
5
```

By providing a `window` and `step` size, counting can be done from within
a sliding window:

```
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

```
julia> count_pairwise(Match, dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA")
4×4 Array{Int64,2}:
 0  6  7  5
 6  0  7  7
 7  7  0  6
 5  7  6  0
```




Recall that `LongDNASeq` is a type alias of `GeneralSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `GeneralSequence{DNAAlphabet{2}}` saves half
of memory footprint compared to `GeneralSequence{DNAAlphabet{4}}`. If you need to
handle reference genomes that are composed of five nucleotides, ACGTN,
consider to use the `ReferenceSequence` type described in the [Reference
sequences](@ref) section.
