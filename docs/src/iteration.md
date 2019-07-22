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
               n += 1
           end
       end

julia> n
3

```

However, there are many other ways to iterate over a biological sequence, and
they are described in the following sections.

## Kmers and Skipmers

To iterate every kmer in a longer DNA or RNA sequence, use the `each` method:

```@docs
each(::Type{Kmer{T,K}}, ::BioSequence, ::Integer)
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

