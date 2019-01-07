```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# General-purpose sequences

`GeneralSequence{A}` is a generic sequence type parameterized by an
[alphabet](@ref alphabet) type `A` that defines the domain (or set) of
biological symbols, and each alphabet has an associated symbol type.
For example, `AminoAcidAlphabet` is associated with `AminoAcid` and hence an
object of the `GeneralSequence{AminoAcidAlphabet}` type represents a sequence of
amino acids.  Symbols from multiple alphabets can't be intermixed in one
sequence type.

The following table summarizes common sequence types that are defined in the
`BioSequences` module:

| Type                                   | Symbol type | Type alias          |
| :------------------------------------- | :---------- | :------------------ |
| `GeneralSequence{DNAAlphabet{4}}`      | `DNA`       | `DNASequence`       |
| `GeneralSequence{RNAAlphabet{4}}`      | `RNA`       | `RNASequence`       |
| `GeneralSequence{AminoAcidAlphabet}`   | `AminoAcid` | `AminoAcidSequence` |
| `GeneralSequence{CharAlphabet}`        | `Char`      | `CharSequence`      |

Parameterized definition of the `GeneralSequence{A}` type is for the purpose of
unifying the data structure and operations of any symbol type. In most cases,
users don't have to care about it and can use *type aliases* listed above.
However, the alphabet type fixes the internal memory encoding and plays an
important role when optimizing performance of a program
(see [Using a more compact sequence representation](@ref) section for low-memory
encodings).  It also enables a developer to define their own alphabet only by
defining few numbers of methods.
This is described in [Defining a new alphabet](@ref) section.


## Constructing `GeneralSequence`s

### Using string literals

Most immediately, sequence literals can be constructed using the string macros
`dna`, `rna`, `aa`, and `char`:

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

julia> char"αβγδϵ"
5char Char Sequence:
αβγδϵ

```

However it should be noted that by default these sequence literals
allocate the `GeneralSequence` object before the code containing the sequence
literal is run.
This means there may be occasions where your program does not behave as you
first expect.
For example consider the following code:

```jldoctest
julia> function foo()
           s = dna"CTT"
           push!(s, DNA_A)
       end
foo (generic function with 1 method)

```

```@meta
DocTestSetup = quote
    using BioSequences
    function foo()
        s = dna"CTT"d
        push!(s, DNA_A)
    end
end
```

You might expect that every time you call `foo`, that a DNA sequence `CTTA` would
be returned. You might expect that this is because every time `foo` is called,
a new DNA sequence variable `CTT` is created, and and `A` nucleotide is pushed
to it, and the result, `CTTA` is returned.
In other words you might expect the following output:

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

```

However, this is not what happens, instead the following happens:

```@meta
DocTestSetup = quote
    using BioSequences
    function foo()
        s = dna"CTT"s
        push!(s, DNA_A)
    end
end
```

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
5nt DNA Sequence:
CTTAA

julia> foo()
6nt DNA Sequence:
CTTAAA

```

The reason for this is because the sequence literal is allocated only once
before the first time the function `foo` is called and run. Therefore, `s` in
`foo` is always a reference to that one sequence that was allocated.
So one sequence is created before `foo` is called, and then it is pushed to
every time `foo` is called. Thus, that one allocated sequence grows with every
call of `foo`.

If you wanted `foo` to create a new sequence each time it is called,
then you can add a flag to the end of the sequence literal to dictate behaviour:
A flag of 's' means 'static': the sequence will be allocated before code is run,
as is the default behaviour described above.
However providing 'd' flag changes the behaviour: 'd' means 'dynamic':
the sequence will be allocated at whilst the code is running, and not before.
So to change `foo` so as it creates a new sequence
each time it is called, simply add the 'd' flag to the sequence literal:
```@meta
DocTestSetup = quote
    using BioSequences
end
```

```jldoctest
julia> function foo()
           s = dna"CTT"d     # 'd' flag appended to the string literal.
           push!(s, DNA_A)
       end
foo (generic function with 1 method)

```

Now every time `foo` is called, a new sequence `CTT` is created, and an `A`
nucleotide is pushed to it:

```@meta
DocTestSetup = quote
    using BioSequences
    function foo()
        s = dna"CTT"d
        push!(s, DNA_A)
    end
end
```

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

```

So the take home message of sequence literals is this:

Be careful when you are using sequence literals inside of functions, and inside
the bodies of things like for loops. And if you use them and are unsure, use the
 's' and 'd' flags to ensure the behaviour you get is the behaviour you intend.


### Other constructors and conversion

`GeneralSequence`s can also be constructed from strings or arrays of nucleotide
or amino acid symbols using constructors or the `convert` function:

```jldoctest
julia> DNASequence("TTANC")
5nt DNA Sequence:
TTANC

julia> DNASequence([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

julia> convert(DNASequence, [DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

```

Using `convert`, these operations are reversible: sequences can be converted to
strings or arrays:
```jldoctest
julia> convert(String, dna"TTANGTA")
"TTANGTA"

julia> convert(Vector{DNA}, dna"TTANGTA")
7-element Array{BioSymbols.DNA,1}:
 DNA_T
 DNA_T
 DNA_A
 DNA_N
 DNA_G
 DNA_T
 DNA_A

```

Sequences can also be concatenated into longer sequences:
```jldoctest
julia> DNASequence(dna"ACGT", dna"NNNN", dna"TGCA")
12nt DNA Sequence:
ACGTNNNNTGCA

julia> dna"ACGT" * dna"TGCA"
8nt DNA Sequence:
ACGTTGCA

julia> repeat(dna"TA", 10)
20nt DNA Sequence:
TATATATATATATATATATA

julia> dna"TA" ^ 10
20nt DNA Sequence:
TATATATATATATATATATA

```

Despite being separate types, `DNASequence` and `RNASequence` can freely be
converted between efficiently without copying the underlying data:
```jldoctest
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> rna = convert(RNASequence, dna)
12nt RNA Sequence:
UUANGUAGACCG

julia> dna.data === rna.data  # underlying data are same
true

```

A random sequence can be obtained by the `randdnaseq`, `randrnaseq` and
`randaaseq` functions, which generate `DNASequence`, `RNASequence` and
`AminoAcidSequence`, respectively. Generated sequences are composed of the
standard symbols without ambiguity and gap. For example, `randdnaseq(6)` may
generate `dna"TCATAG"` but never generates `dna"TNANAG"` or `dna"T-ATAG"`.

A translatable `RNASequence` can also be converted to an `AminoAcidSequence`
using the [`translate`](@ref) function.






### Setindex and modifying DNA sequences





## Site counting

BioSequences extends the `Base.count` method to provide some useful utilities for
counting the number of sites in biological sequences.

### Site types

Different types of site can be counted. Each of these types is a concrete
subtype of the abstract type `Site`:

```@docs
Certain
Gap
Ambiguous
Match
Mismatch
```

### `Base.count` methods

The count method can be used with two sequences and a concrete subtype of
`Site`:

```jldoctest
julia> count(Match, dna"ATCGATCG", dna"AAGGTTCG")
5
```

By providing a `window` and `step` size, counting can be done from within
a sliding window:

```jldoctest
julia> count(Match, dna"ATCGATCG", dna"AAGGTTCG", 3, 1)
6-element Array{IntervalTrees.IntervalValue{Int64,Int64},1}:
 IntervalTrees.IntervalValue{Int64,Int64}
(1,3) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(2,4) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(3,5) => 1
 IntervalTrees.IntervalValue{Int64,Int64}
(4,6) => 2
 IntervalTrees.IntervalValue{Int64,Int64}
(5,7) => 2
 IntervalTrees.IntervalValue{Int64,Int64}
(6,8) => 3
```

### The `pairwise_count` function
Counting can also be done on a set of sequences in a pairwise manner with the
`count_pairwise` function:

```jldoctest
julia> count_pairwise(Match, dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA")
4×4 Array{Int64,2}:
 0  6  7  5
 6  0  7  7
 7  7  0  6
 5  7  6  0
```

## Iteration

Sequences also work as iterators over symbols:

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


## Using a more compact sequence representation

As we saw above, DNA and RNA sequences can store any ambiguous nucleotides like
'N'.
If you are sure that nucleotide sequences store unambiguous nucleotides
only, you can save the memory space of sequences by using a slightly different
type:
`DNAAlphabet{2}` is an alphabet that uses two bits per base and limits to only
unambiguous nucleotide symbols (ACGT in DNA and ACGU in RNA).
To create a sequence of this alphabet, you need to explicitly pass
`DNAAlphabet{2}` to `BioSequence` as its type parameter:

```jldoctest
julia> seq = GeneralSequence{DNAAlphabet{2}}("ACGT")
4nt DNA Sequence:
ACGT

```

Recall that `DNASequence` is a type alias of `GeneralSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `GeneralSequence{DNAAlphabet{2}}` saves half
of memory footprint compared to `GeneralSequence{DNAAlphabet{4}}`. If you need to
handle reference genomes that are composed of five nucleotides, ACGTN,
consider to use the `ReferenceSequence` type described in the [Reference
sequences](@ref) section.
