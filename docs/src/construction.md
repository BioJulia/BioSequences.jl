```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Construction & conversion

Here we will showcase the various ways you can construct the various sequence
types in BioSequences.

## Constructing sequences

### From strings

Sequences can be constructed from strings using their constructors:

```jldoctest
julia> LongDNA{4}("TTANC")
5nt DNA Sequence:
TTANC

julia> LongSequence{DNAAlphabet{2}}("TTAGC")
5nt DNA Sequence:
TTAGC

julia> LongRNA{4}("UUANC")
5nt RNA Sequence:
UUANC

julia> LongSequence{RNAAlphabet{2}}("UUAGC")
5nt RNA Sequence:
UUAGC

```

Type alias' can also be used for brevity.

```jldoctest
julia> LongDNA{4}("TTANC")
5nt DNA Sequence:
TTANC

julia> LongDNA{2}("TTAGC")
5nt DNA Sequence:
TTAGC

julia> LongRNA{4}("UUANC")
5nt RNA Sequence:
UUANC

julia> LongRNA{2}("UUAGC")
5nt RNA Sequence:
UUAGC
```

### Constructing sequences from arrays of BioSymbols

Sequences can be constructed using vectors or arrays of a `BioSymbol` type:

```jldoctest
julia> LongDNA{4}([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

julia> LongSequence{DNAAlphabet{2}}([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
5nt DNA Sequence:
TTAGC

```

### Constructing sequences from other sequences

You can create sequences, by concatenating other sequences together:
```jldoctest
julia> LongDNA{2}("ACGT") * LongDNA{2}("TGCA")
8nt DNA Sequence:
ACGTTGCA

julia> repeat(LongDNA{4}("TA"), 10)
20nt DNA Sequence:
TATATATATATATATATATA

julia> LongDNA{4}("TA") ^ 10
20nt DNA Sequence:
TATATATATATATATATATA

```

Sequence views (`LongSubSeq`s) are special, in that they do not own their own data,
and must be constructed from a `LongSequence` or another `LongSubSeq`:

```jdoctest
julia> seq = LongDNA{4}("TACGGACATTA")
11nt DNA Sequence:
TACGGACATTA

julia> seqview = LongSubSeq(seq, 3:7)
5nt DNA Sequence:
CGGAC

julia> seqview2 = @view seq[1:3]
3nt DNA Sequence:
TAC

julia> typeof(seqview) == typeof(seqviev2) && typeof(seqview) <: LongSubSeq
true

```


## Conversion of sequence types

You can convert between sequence types.
```jldoctest
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> rna = convert(LongRNA{4}, dna)
12nt RNA Sequence:
UUANGUAGACCG

```

Sequences can be converted explicitly and implicitly, into arrays and strings:

```jldoctest
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> dnastr = convert(String, dna)
"TTANGTAGACCG"

julia> # Implicit conversion to string - putting dna sequence in String vector 

julia> arr = String[dna]
1-element Vector{String}:
 "TTANGTAGACCG"
 
```

## String literals

BioSequences provides several string literal macros for creating sequences.

!!! note
    When you use literals you may mix the case of characters.

### Long sequence literals

```jldoctest
julia> dna"TACGTANNATC"
11nt DNA Sequence:
TACGTANNATC

julia> rna"AUUUGNCCANU"
11nt RNA Sequence:
AUUUGNCCANU

julia> aa"ARNDCQEGHILKMFPSTWYVX"
21aa Amino Acid Sequence:
ARNDCQEGHILKMFPSTWYVX
```

Unlike most string literal macros in Julia, the resulting sequences are allocated
at runtime, not at compile time - so called _dynamic allocation_. Most other Julia
string macros create the objects at compile time.

You can create biosequences at compile time by setting the "static" flag `s` after
the literal. This may improve performance by moving the allocation to compile time,
but can lead to unexpected behaviour, as the following example shows:

```jldoctest
julia> f() = dna"C"s;

julia> push!(f(), DNA_A)
2nt DNA Sequence:
CA

julia> push!(f(), DNA_A) # same sequence - two A's!
3nt DNA Sequence:
CAA
```

In the example above, the same `const` (static) sequence is being returned from `f()`,
and grown at each `push!` call. This makes `f` a fast function, since it does not need
to allocate memory.

The default behaviour, namely dynamic allocation at runtime, can also be obtained by
using the `d` flag:

```jldoctest
julia> f() = dna"C"d;

julia> push!(f(), DNA_A)
2nt DNA Sequence:
CA

julia> push!(f(), DNA_A) # new sequence - one A!
2nt DNA Sequence:
CA
```
