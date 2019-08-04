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

You can also use two sequences, for example to compute the number of matching
or mismatching symbols:

```jldoctest
julia> count(!=, dna"ATCGM", dna"GCCGM")
2

julia> count(==, dna"ATCGM", dna"GCCGM")
3

```

## Alias functions

A number of functions which are aliases for various invocations of `Base.count`
are provided.

| Alias function | Base.count call(s)                                          |
| :------------- | :---------------------------------------------------------- |
| `n_ambiguous`  | `count(isambiguous, seq)`, `count(isambiguous, seqa, seqb)` |
| `n_certain`    | `count(iscertain, seq)`, `count(iscertain, seqa, seqb)`     |
| `n_gap`        | `count(isgap, seq)`, `count(isgap, seqa, seqb)`             |
| `matches`      | `count(==, seqa, seqb)`                                     |
| `mismatches`   | `count(!=, seqa, seqb)`                                     |

## Bit-parallel optimisations

For the vast majority of `Base.count(f, seq)` and `Base.count(f, seqa, seqb)`
methods, a naive counting is done: the internal `count_naive` function is called,
which simply loops over each position, applies `f`, and accumulates the result.

However, for some functions, it is possible to implement highly efficient methods
that use bit-parallelism to check many elements at one time.
This is made possible by the succinct encoding of BioSequences.
Usually `f` is one of the functions provided by BioSymbols.jl or by BioSequences.jl

For such sequence and function combinations, `Base.count(f, seq)` is overloaded
to call an internal `BioSequences.count_*_bitpar` function, which is passed the
sequence(s). If you want to force BioSequences to use naive counting for the
purposes of testing or debugging for example, then you can call
`BioSequences.count_naive` directly.