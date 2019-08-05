```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Indexing & modifying sequences

## Indexing

Most `BioSequence` concrete subtypes for the most part behave like other vector
or string types. They can be indexed using integers or ranges:

For example, with `LongSequence`s:

```jldoctest
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5]
DNA_T

julia> seq[6:end]
14nt DNA Sequence:
TANAGTNNAGTACC

```

The biological symbol at a given locus in a biological sequence can be set using
setindex:

```jldoctest
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5] = DNA_A
DNA_A

```

!!! note
    Some types such as `Kmer` can be indexed using integers but not using ranges.

For `LongSequence` types, indexing a sequence by range creates a subsequence
of the original sequence.
Unlike `Arrays` in the standard library, creating this subsequence is copy-free:
the subsequence simply points to the original `LongSequence` data with its range.
You may think that this is unsafe because modifying subsequences propagates to
the original sequence, but this doesn't happen actually:

```jldoctest
julia> seq = dna"AAAA"    # create a sequence
4nt DNA Sequence:
AAAA

julia> subseq = seq[1:2]  # create a subsequence from `seq`
2nt DNA Sequence:
AA

julia> subseq[2] = DNA_T  # modify the second element of it
DNA_T

julia> subseq             # the subsequence is modified
2nt DNA Sequence:
AT

julia> seq                # but the original sequence is not
4nt DNA Sequence:
AAAA

```

This is because modifying a sequence checks whether its underlying data are
shared with other sequences under the hood.
If and only if the data are shared, the subsequence creates a copy of itself.
Any modifying operation does this check.
This is called *copy-on-write* strategy and users don't need to care
about it because it is transparent: If the user modifies a sequence with or
subsequence, the job of managing and protecting the underlying data of sequences
is handled for them.

## Modifying sequences

In addition to `setindex`, many other modifying operations are possible for
biological sequences such as `push!`, `pop!`, and `insert!`, which should be
familiar to anyone used to editing arrays.

```@docs
push!
pop!
pushfirst!
popfirst!
insert!
deleteat!(::BioSequences.BioSequence, ::Integer)
append!
resize!
empty!
```

Here are some examples:

```jldoctest
julia> seq = dna"ACG"
3nt DNA Sequence:
ACG

julia> push!(seq, DNA_T)
4nt DNA Sequence:
ACGT

julia> append!(seq, dna"AT")
6nt DNA Sequence:
ACGTAT

julia> deleteat!(seq, 2)
5nt DNA Sequence:
AGTAT

julia> deleteat!(seq, 2:3)
3nt DNA Sequence:
AAT

```

### Additional transformations

In addition to these basic modifying functions, other sequence transformations
that are common in bioinformatics are also provided.

```@docs
reverse!(::BioSequences.LongSequence)
reverse(::BioSequences.LongSequence{<:NucleicAcidAlphabet})
complement!
complement
reverse_complement!
reverse_complement
ungap!
ungap
canonical!
canonical
```

Some examples:

```jldoctest
julia> seq = dna"ACGTAT"
6nt DNA Sequence:
ACGTAT

julia> reverse!(seq)
6nt DNA Sequence:
TATGCA

julia> complement!(seq)
6nt DNA Sequence:
ATACGT

julia> reverse_complement!(seq)
6nt DNA Sequence:
ACGTAT

```

Many of these methods also have a version which makes a copy of the input
sequence, so you get a modified copy, and don't alter the original sequence.
Such methods are named the same, but without the exclamation mark.
E.g. `reverse` instead of `reverse!`, and `ungap` instead of `ungap!`.  

#### Translation

Translation is a slightly more complex transformation for RNA Sequences and so
we describe it here in more detail.

The [`translate`](@ref) function translates a sequence of codons in a RNA sequence
to a amino acid sequence based on a genetic code. The `BioSequences` package
provides all NCBI defined genetic codes and they are registered in
[`ncbi_trans_table`](@ref).

```@docs
translate
ncbi_trans_table
```

```jldoctest
julia> ncbi_trans_table
Translation Tables:
  1. The Standard Code (standard_genetic_code)
  2. The Vertebrate Mitochondrial Code (vertebrate_mitochondrial_genetic_code)
  3. The Yeast Mitochondrial Code (yeast_mitochondrial_genetic_code)
  4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (mold_mitochondrial_genetic_code)
  5. The Invertebrate Mitochondrial Code (invertebrate_mitochondrial_genetic_code)
  6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (ciliate_nuclear_genetic_code)
  9. The Echinoderm and Flatworm Mitochondrial Code (echinoderm_mitochondrial_genetic_code)
 10. The Euplotid Nuclear Code (euplotid_nuclear_genetic_code)
 11. The Bacterial, Archaeal and Plant Plastid Code (bacterial_plastid_genetic_code)
 12. The Alternative Yeast Nuclear Code (alternative_yeast_nuclear_genetic_code)
 13. The Ascidian Mitochondrial Code (ascidian_mitochondrial_genetic_code)
 14. The Alternative Flatworm Mitochondrial Code (alternative_flatworm_mitochondrial_genetic_code)
 16. Chlorophycean Mitochondrial Code (chlorophycean_mitochondrial_genetic_code)
 21. Trematode Mitochondrial Code (trematode_mitochondrial_genetic_code)
 22. Scenedesmus obliquus Mitochondrial Code (scenedesmus_obliquus_mitochondrial_genetic_code)
 23. Thraustochytrium Mitochondrial Code (thraustochytrium_mitochondrial_genetic_code)
 24. Pterobranchia Mitochondrial Code (pterobrachia_mitochondrial_genetic_code)
 25. Candidate Division SR1 and Gracilibacteria Code (candidate_division_sr1_genetic_code)

```

<https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes>