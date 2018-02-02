```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Sequence search

Three kinds of on-line search functions are provided:

1. Exact search
2. Approximate search
3. Regular expression search

These are all specialized for biological sequences and ambiguities of symbols
are considered.

## Exact search

Exact search functions search for an occurrence of the query symbol or
sequence.
```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> findfirst(DNA_G, seq)
4

julia> query = dna"AGC";

julia> findfirst(query, seq)
3:5

julia> findlast(query, seq)
8:10

```

These search functions take ambiguous symbols into account.
That is, if two symbols are compatible (e.g. `DNA_A` and `DNA_N`),
they match when searching an occurrence.
In the following example, 'N' is a wild card that matches any symbols.

```jldoctest
julia> findfirst(dna"CGT", dna"ACNT")  # 'N' matches 'G'
2:4

julia> findfirst(dna"CNT", dna"ACGT")  # 'G' matches 'N'
2:4

```

The exception to this behaviour is if you are finding a single 'character',
in which case an ambiguous symbol is matched exactly:

```jldoctest
julia> findfirst(DNA_N, dna"ACNT")
3

```

The exact sequence search needs a preprocessing phase of query sequence before
the searching phase. This would be fast enough for most search applications.
But when searching a query sequence to many target sequences, caching
the result of preprocessing may save time. You can do this by creating an
`ExactSearchQuery` object and re-use it for each search:
```jldoctest
julia> query = ExactSearchQuery(dna"ATT");

julia> findfirst(query, dna"ATTTATT")
1:3

julia> findlast(query, dna"ATTTATT")
5:7

```


## Approximate search

The approximate search is similar to the exact search but allows a specific
number of errors. That is, it tries to find a subsequence of the target sequence
within a specific [Levenshtein
distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of the query
sequence:
```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> approxsearch(seq, dna"AGGG", 0)  # nothing matches with no errors
0:-1

julia> approxsearch(seq, dna"AGGG", 1)  # seq[3:5] matches with one error
3:6

julia> approxsearch(seq, dna"AGGG", 2)  # seq[1:4] matches with two errors
1:4

```

Like the exact search functions, four kinds of functions (`approxsearch`,
`approxsearchindex`, `approxrsearch`, and `approxrsearchindex`) are available:
```jldoctest
julia> seq = dna"ACAGCGTAGCT"; pat = dna"AGGG";

julia> approxsearch(seq, pat, 2)        # return the range (forward)
1:4

julia> approxsearchindex(seq, pat, 2)   # return the starting index (forward)
1

julia> approxrsearch(seq, pat, 2)       # return the range (backward)
8:11

julia> approxrsearchindex(seq, pat, 2)  # return the starting index (backward)
8

```

Preprocessing can be cached in an `ApproximateSearchQuery` object:
```jldoctest
julia> query = ApproximateSearchQuery(dna"AGGG");

julia> approxsearch(dna"AAGAGG", query, 1)
2:5

julia> approxsearch(dna"ACTACGT", query, 2)
4:6

```

## Regular expression search

Query patterns can be described in regular expressions. The syntax supports
a subset of Perl and PROSITE's notation.

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

`match` always returns a `Nullable` object and it should be null if no match is
found.

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

`eachmatch` and `findfirst` are also defined like usual strings:

```jldoctest
julia> collect(matched(x) for x in eachmatch(biore"TATA*?"d, dna"TATTATAATTA")) # overlap
4-element Array{BioSequence{DNAAlphabet{4}},1}:
 TAT  
 TAT
 TATA
 TATAA

julia> collect(matched(x) for x in eachmatch(biore"TATA*"d, dna"TATTATAATTA", false)) # no overlap
2-element Array{BioSequence{DNAAlphabet{4}},1}:
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


## Position weight matrix search

A motif can also be specified using [position weight
matrix](https://en.wikipedia.org/wiki/Position_weight_matrix) (PWM) in a
probabilistic way. `search(seq, pwm, threshold)` method searches for the first
position in the sequence where a score calculated using the PWM is greater than
or equal to the threshold. More formally, denoting the sequence as ``S`` and the
PWM value of symbol ``s`` at position ``j`` as ``M_{s,j}``, the score starting
from a position ``p`` is defined as

```math
\operatorname{score}(S, p) = \sum_{i=1}^L M_{S[p+i-1],i}
```

and `search(S, M, t)` returns the smallest ``p`` that satisfies
``\operatorname{score}(S, p) \ge t``.

There are two kinds of matrices in this package: `PFM` and `PWM`. The `PFM` type
is a position frequency matrix and stores symbol frequencies for each position.
The `PWM` is a position weight matrix and stores symbol scores for each
position. You can create a `PFM` from a set of sequences with the same length
and then create a `PWM` from the `PFM` object.

```jldoctest
julia> kmers = DNAKmer.(["TTA", "CTA", "ACA", "TCA", "GTA"])
5-element Array{Kmer{DNA,3},1}:
 TTA
 CTA
 ACA
 TCA
 GTA

julia> pfm = PFM(kmers)  # sequence set => PFM
4×3 PFM{DNA,Int64}:
 A  1  0  5
 C  1  2  0
 G  1  0  0
 T  2  3  0

julia> pwm = PWM(pfm)  # PFM => PWM
4×3 PWM{DNA,Float64}:
 A -0.321928 -Inf       2.0
 C -0.321928  0.678072 -Inf
 G -0.321928 -Inf      -Inf
 T  0.678072  1.26303  -Inf

julia> pwm = PWM(pfm .+ 0.01)  # add pseudo counts to avoid infinite values
4×3 PWM{DNA,Float64}:
 A -0.319068 -6.97728   1.99139
 C -0.319068  0.673772 -6.97728
 G -0.319068 -6.97728  -6.97728
 T  0.673772  1.25634  -6.97728

julia> pwm = PWM(pfm .+ 0.01, prior=[0.2, 0.3, 0.3, 0.2])  # GC-rich prior
4×3 PWM{DNA,Float64}:
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

[Wasserman2004]: https://doi.org/10.1038/nrg1315
