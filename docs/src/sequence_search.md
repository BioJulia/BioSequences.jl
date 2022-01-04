```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Searching for sequence motifs

There are many ways to search for particular motifs in biological sequences:

1. Exact searches, where you are looking for exact matches of a particular
   character of substring.
2. Approximate searches, where you are looking for sequences that are
   sufficiently similar to a given sequence or family of sequences.
3. Searches where you are looking for sequences that conform to some sort of
   pattern.

All these kinds of searches are provided in BioSequences.jl, and they all 
conform to the `findnext`, `findprev`, and `occursin` patterns established in `Base` for
`String` and collections like `Vector`. The exception is searching using the specialised
regex provided in this package, which as you shall see, conforms to the `match`
pattern established in `Base` for pcre and `String`s.


## Exact search

Exact searches, are where you are looking for exact matches of a particular
character of substring.

Like other Julia sequences such as `Vector`, you can search a `BioSequence` with the `findfirst(predicate, collection)` method pattern. 
The method applies the predicate function to elements in the collection, and when the predicate evaluates as `true`, the method returns the index of the element. 
If the predicate does not evaluate as `true` for any elements in the collection, the method returns `nothing`.

For example, the following will find the first exact match of a specific
symbol.

```jldoctest
julia> findfirst(isequal(DNA_A), dna"GCTTAG")
5

julia> findfirst(isequal(DNA_M), dna"GCTTAG") === nothing
true
```

To search for a subsequence, first, construct an `ExactSearchQuery` to pass to
the findX functions.

```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> query = ExactSearchQuery(dna"AGC");

julia> findfirst(query, seq)
3:5

julia> findlast(query, seq)
8:10

julia> occursin(query, seq)
true
```

By default, `ExactSearchQuery` uses the `isequal` function to compare symbols. 
However, in biology, sometimes we want a more flexible comparison to find subsequences of _compatible_ symbols. 
To do this, we can specify the `iscompatible` function in the constructor of `ExactSearchQuery`.

```jldoctest
julia> findfirst(ExactSearchQuery(dna"CGT", iscompatible), dna"ACNT")  # 'N' matches 'G'
2:4

julia> findfirst(ExactSearchQuery(dna"CNT", iscompatible), dna"ACGT")  # 'G' matches 'N'
2:4

julia> occursin(ExactSearchQuery(dna"CNT", iscompatible), dna"ACNT")
true
```


## Allowing mismatches

The approximate search is similar to the exact search but allows a specific
number of errors. That is, it tries to find a subsequence of the target sequence
within a specific [Levenshtein
distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of the query
sequence:

```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> query = ApproximateSearchQuery(dna"AGGG");

julia> findfirst(query, 0, seq)  # nothing matches with no errors
nothing

julia> findfirst(query, 1, seq)  # seq[3:6] matches with one error
3:6

julia> findfirst(query, 2, seq)  # seq[1:4] matches with two errors
1:4

```

Like the `ExactSearchQuery`, all the findX functions including `findfirst`,
`findlast`, `findnext`, `findprev`, and `occursin` can be used with
`ApproximateSearchQuery`, you simply have to provide the query as the first parameter,
and the number of mismatches you will allow as the second.

```jldoctest
julia> query = ApproximateSearchQuery(dna"AGGG");

julia> occursin(query, 1, dna"AAGAGG")
true

julia> findnext(query, 2, dna"ACTACGT", 1)
4:6

```

Also, like the `ExactSearchQuery`, you can pass a comparator function such as
`isequal` or `iscompatible` to its constructor to modify the search behaviour.

```jldoctest
julia> query = ApproximateSearchQuery(dna"AGGG", iscompatible);

julia> occursin(query, 1, dna"AAGNGG")    # 1 mismatch permitted (A vs G) & matched N
true

julia> findnext(query, 1, dna"AAGNGG", 1) # 1 mismatch permitted (A vs G) & matched N
1:4

```

!!! note
   This method of searching for motifs was implemented with smaller query motifs
   in mind.
   
   If you are looking to search for imperfect matches of longer sequences in this
   manner, you are likely better off using some kind of local-alignment algorithm
   or one of the BLAST variants.

## Searching according to a pattern

### Regular expression search

Query patterns can be described in regular expressions. The syntax supports
a subset of Perl and PROSITE's notation.

Biological regexes can be constructed using the `BioRegex` constructor, for
example by doing `BioRegex{AminoAcid}("MV+")`. For bioregex literals, it is
instead recommended using the `@biore_str` macro:

The Perl-like syntax starts with `biore` (BIOlogical REgular expression)
and ends with a symbol option: "dna", "rna" or "aa". For example, `biore"A+"dna`
is a regular expression for DNA sequences and `biore"A+"aa` is for amino acid
sequences. The symbol options can be abbreviated to its first character: "d",
"r" or "a", respectively.

Here are examples of using the regular expression for `BioSequence`s:
```jldoctest
julia> match(biore"A+C*"dna, dna"AAAACC")
RegexMatch("AAAACC")

julia> match(biore"A+C*"d, dna"AAAACC")
RegexMatch("AAAACC")

julia> occursin(biore"A+C*"dna, dna"AAC")
true

julia> occursin(biore"A+C*"dna, dna"C")
false

```

`match` will return a `RegexMatch` if a match is found, otherwise it will return `nothing` if no match is found.

The table below summarizes available syntax elements.

| Syntax | Description | Example |
|:------:|:------------|:--------|
| `\|` | alternation | `"A\|T"` matches `"A"` and `"T"` |
| `*` | zero or more times repeat | `"TA*"` matches `"T"`, `"TA"` and `"TAA"` |
| `+` | one or more times repeat | `"TA+"` matches `"TA"` and `"TAA"` |
| `?` | zero or one time | `"TA?"` matches `"T"` and `"TA"` |
| `{n,}` | `n` or more times repeat | `"A{3,}"` matches `"AAA"` and `"AAAA"` |
| `{n,m}` | `n`-`m` times repeat | `"A{3,5}"` matches `"AAA"`, `"AAAA"` and `"AAAAA"`|
| `^` | the start of the sequence | `"^TAN*"` matches `"TATGT"` |
| `$` | the end of the sequence | `"N*TA$"` matches `"GCTA"` |
| `(...)` | pattern grouping | `"(TA)+"` matches `"TA"` and `"TATA"` |
| `[...]` | one of symbols | `"[ACG]+"` matches `"AGGC"` |

`eachmatch` and `findfirst` are also defined, just like usual regex and strings
found in `Base`.

```jldoctest
julia> collect(matched(x) for x in eachmatch(biore"TATA*?"d, dna"TATTATAATTA")) # overlap
4-element Vector{LongSequence{DNAAlphabet{4}}}:
 TAT  
 TAT
 TATA
 TATAA

julia> collect(matched(x) for x in eachmatch(biore"TATA*"d, dna"TATTATAATTA", false)) # no overlap
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 TAT  
 TATAA

julia> findfirst(biore"TATA*"d, dna"TATTATAATTA")
1:3

julia> findfirst(biore"TATA*"d, dna"TATTATAATTA", 2)
4:8

```

Noteworthy differences from strings are:

* Ambiguous characters match any compatible characters (e.g. `biore"N"d` is equivalent to `biore"[ACGT]"d`).
* Whitespaces are ignored (e.g. `biore"A C G"d` is equivalent to `biore"ACG"d`).

The PROSITE notation is described in [ScanProsite - user
manual](https://prosite.expasy.org/scanprosite/scanprosite_doc.html). The syntax
supports almost all notations including the extended syntax. The PROSITE
notation starts with `prosite` prefix and no symbol option is needed because it
always describes patterns of amino acid sequences:

```jldoctest
julia> match(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVPQARG")
RegexMatch("CPVPQARG")

julia> match(prosite"[AC]xVx(4){ED}", aa"CPVPQARG")
RegexMatch("CPVPQARG")

```


### Position weight matrix search

A motif can be specified using [position weight
matrix](https://en.wikipedia.org/wiki/Position_weight_matrix) (PWM) in a
probabilistic way. 
This method searches for the first position in the sequence where a score
calculated using a PWM is greater than or equal to a threshold.
More formally, denoting the sequence as ``S`` and the PWM value of symbol ``s``
at position ``j`` as ``M_{s,j}``, the score starting from a position ``p`` is
defined as

```math
\operatorname{score}(S, p) = \sum_{i=1}^L M_{S[p+i-1],i}
```

and the search returns the smallest ``p`` that satisfies
``\operatorname{score}(S, p) \ge t``.

There are two kinds of matrices in this package: `PFM` and `PWM`. The `PFM` type
is a position frequency matrix and stores symbol frequencies for each position.
The `PWM` is a position weight matrix and stores symbol scores for each
position. You can create a `PFM` from a set of sequences with the same length
and then create a `PWM` from the `PFM` object.

```jldoctest
julia> motifs = [dna"TTA", dna"CTA", dna"ACA", dna"TCA", dna"GTA"]
5-element Vector{LongSequence{DNAAlphabet{4}}}:
 TTA
 CTA
 ACA
 TCA
 GTA

julia> pfm = PFM(motifs)  # sequence set => PFM
4×3 PFM{DNA, Int64}:
 A  1  0  5
 C  1  2  0
 G  1  0  0
 T  2  3  0

julia> pwm = PWM(pfm)  # PFM => PWM
4×3 PWM{DNA, Float64}:
 A -0.321928 -Inf       2.0
 C -0.321928  0.678072 -Inf
 G -0.321928 -Inf      -Inf
 T  0.678072  1.26303  -Inf

julia> pwm = PWM(pfm .+ 0.01)  # add pseudo counts to avoid infinite values
4×3 PWM{DNA, Float64}:
 A -0.319068 -6.97728   1.99139
 C -0.319068  0.673772 -6.97728
 G -0.319068 -6.97728  -6.97728
 T  0.673772  1.25634  -6.97728

julia> pwm = PWM(pfm .+ 0.01, prior=[0.2, 0.3, 0.3, 0.2])  # GC-rich prior
4×3 PWM{DNA, Float64}:
 A  0.00285965 -6.65535   2.31331
 C -0.582103    0.410737 -7.24031
 G -0.582103   -7.24031  -7.24031
 T  0.9957      1.57827  -6.65535

```

The ``PWM_{s,j}`` matrix is computed from ``PFM_{s,j}`` and the prior
probability ``p(s)`` as follows ([Wasserman2004]):

```math
\begin{align}
    PWM_{s,j} &= \log_2 \frac{p(s,j)}{p(s)} \\
    p(s,j)  &= \frac{PFM_{s,j}}{\sum_{s'} PFM_{s',j}}.
\end{align}
```

However, if you just want to quickly conduct a search, constructing the pfm and
pwm is done for you as a convenience if you build a `PWMSearchQuery`, using a
collection of sequences:

```jldoctest
subject = dna"TATTATAATTA"
qa = PWMSearchQuery(motifs)
findfirst(qa, subject)
```

[Wasserman2004]: https://doi.org/10.1038/nrg1315
