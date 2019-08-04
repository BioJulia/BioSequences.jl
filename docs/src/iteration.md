```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Iteration

As you might expect, sequence types are iterators over their elements:

```jldoctest
julia> n = 0
0

julia> for nt in dna"ATNGNNT"
           if nt == DNA_N
               global n += 1
           end
       end

julia> n
3

```

However, there are many other ways to iterate over a biological sequence, and
they are described in the following sections.

## Kmers and Skipmers

### Kmers

To iterate over every overlapping kmer in a longer DNA or RNA sequence, use the
`each` method:

```@docs
each(::Type{T}, seq::BioSequence) where {T<:AbstractMer}
```

Each iteration yields a `MerIterResult` struct that has the following fields:

- `position`: the position in the sequence at which the mer began.
- `fw`: the mer in the same orientation as the sequence from which it was generated.
- `bw`: the reverse complement of the `fw` mer.

Iterating over mers in a sequence builds both a mer and it's reverse complement
at the same time, as it is more efficient, and it is a common requirement; for
example during mapping or constructing assemblies it is often useful to know
if the mer your are processing is a canonical mer or not.

### Kmers with jumps

You can also iterate over non-overlapping kmers using a `step` parameter.

```jlcon
each(DNAMer{27}, seq, 10)
```

### Skipmers

To iterate over Mers using the Skipmer method of selecting nucleotides, then
you provide a `bases_per_cycle` and a `cycle_len` parameter. For example, where
`bases_per_cycle = 2` and `cycle_len = 3`.

```jlcon
each(DNAMer{27}, seq, 2, 3)
```

## Positions

You can iterate over the positions/residues of a sequence that satisfy some
condition or query function:

```@docs
each(::Function, ::BioSequence)
```

```jldoctest
julia> dna_seq = dna"NATTCGRATY"
10nt DNA Sequence:
NATTCGRATY

julia> # Iterate over each ambiguous residue

julia> for (pos, nuc) in each(isambiguous, dna_seq)
           println("Position $pos is $nuc")
       end
Position 1 is N
Position 7 is R
Position 10 is Y

julia> # Iterate over each non-ambiguous residue

julia> for (pos, nuc) in each(iscertain, dna_seq)
           println("Position $pos is $nuc")
       end
Position 2 is A
Position 3 is T
Position 4 is T
Position 5 is C
Position 6 is G
Position 8 is A
Position 9 is T

```

