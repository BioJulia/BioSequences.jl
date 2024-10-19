```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Counting
`BioSequences` contain functionality to efficiently count biosymbols in a biosequence
that satisfies some predicate.

Consider a naive counting function like this:
```julia
function count_Ns(seq::BioSequence{<:DNAAlphabet})
    ns = 0
    for i in seq
        ns += (i == DNA_N)::Bool
    end
    ns
end 
```

This function can be more efficient implemented by exploiting the internal
data encoding of certain biosequences.
Therefore, Julia provides optimised methods for `Base.count`, such that `count_Ns`
above can be more efficiently expressed `count(==(DNA_N), seq)`.

!!! note
    It is important to understand that this speed is achieved with custom methods of
    `Base.count`, and not by a generic mechanism that improves the speed of counting
    symbols in `BioSequence `in general.
    Hence, while `count(==(DNA_N), seq)` may be optimised,
    `count(i -> i == DNA_N, seq)` is not, as this is a different method.

## Currently optimised methods
By default, only `LongSequence` and `LongSubSeq` have optimised methods. Downstream
implementers of new `BioSequence`s may not be optimised.

* `count(isGC, seq)`
* `count(isambiguous, seq)`
* `count(iscertain, seq)`
* `count(isgap, seq)`
* `count(==(biosymbol), seq)` and `count(isequal(biosymbol), seq)`

## Matches and mismatches
The methods `matches` and `mismatches` take two
sequences and count the number of positions where the sequences are unequal or equal, respectively.

They are equivalent to `matches(a, b) = count(splat(==), zip(a, b))`
(and with `!=`, respectively).

```@docs
matches
mismatches
```

## GC content
The convenience function `gc_content(seq)` is equivalent to `count(isGC, seq) / length(seq)`:

```@docs
gc_content
```

## Deprecated aliases
Several of the optimised `count` methods have function names, which are deprecated:

| Deprecated function  | Instead use               |
| :------------------- | :------------------------ |
| `n_gaps`             | `count(isgap, seq)`       |
| `n_certain`          | `count(iscertain, seq)`   |
| `n_ambiguous`        | `count(isambiguous, seq)` |


```@docs
n_gaps
n_certain
n_ambiguous
```
