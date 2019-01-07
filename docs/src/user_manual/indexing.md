```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Indexing

## getindex

Most `BioSequence` concrete subtypes for the most part behave like other vector
or string types. They can be indexed using integers or ranges:

For example, with `GeneralSequence`s:

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

*Note that some types such as `Kmer` can be indexed using integers but not using
ranges.*

For most BioSequence types, indexing a sequence by range creates a subsequence
of the original sequence.
Unlike `Arrays` in the standard library, creating a subsequence is copy-free:
a subsequence simply points to the original sequence
data with its range. You may think that this is unsafe because modifying
subsequences propagates to the original sequence, but this doesn't happen
actually:

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
is handled for them. More details of the *copy-on-write* mechanism are provided
in the developer documentation.

## setindex

The biological symbol at a given locus in a biological sequence can be set using
setindex:

```jldoctest
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5] = DNA_A
DNA_A

```
