```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Sequence composition

There are many instances in analyzing sequence data where you will want to
know about the composition of your sequences.

For example, for a given sequence, you may want to count how many of each
possible Kmer, is present in the sequence. This would be important if - for
instance - you wanted to analyze the Kmer spectra of your data.
Alternatively you might have a collection of sequences, and may want to count
how many of each unique sequence you have in your collection. This would be
important if - for instance - your collection of sequences were from a
population sample, and you wanted to compute the allele or genotype frequencies
for the population.

Whatever the application, BioSequences provides a method called `composition`,
and a parametric struct called `Composition` to both compute, and handle the
results of such sequence composition calculations.

```@docs
Composition{T}
composition
```

Composition structs behave like an associative collection, such as a `Dict`.
But there are a few differences:

1. The `getindex` method for Composition structs is overloaded to return a default
   value of 0, if a key is used that is not present in the Composition.
2. The `merge` method for two Composition structs adds counts together, unlike
   the `merge` method for other associative containers, which would overwrite
   the counts.
