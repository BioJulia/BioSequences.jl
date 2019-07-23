```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Counting

BioSequences extends the `Base.count` method to provide some useful utilities for
counting the number of sites in biological sequences.

Most generically you can count the number of sites that satisfy some condition
i.e. cause some function to return `true`:

```jldoctest
julia> count(isambiguous, dna"ATCGM")
1

```

## Bit-parallel operations

For the vast majority of `count(f, seq)` methods, a naive counting is done:
the internal `count_naive` function is called, which simply loops over each
element of the sequence, applies `f`, and accumulates the result.

However, for some combinations of `f` and `seq`, it is possible to implement
highly efficient methods that use bit-parallelism to check many elements at
one time, which is made possible by the succinct encoding of BioSequences.

Some of these methods have been specialised to take advantage of bit-parallel
operations that process many nucleotides at once.

For such sequence and function combinations, `count(f, seq)` is overloaded to
call the internal `count_bitpar` function, which is passed the sequence, and 
a bit-parallel version of `f`. If you want to force BioSequences to use naive
counting for the purposes of testing or debugging for example, then you can
call `BioSequences.count_naive` directly. 