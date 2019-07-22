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
julia> LongDNASeq("TTANC")
5nt DNA Sequence:
TTANC

julia> LongSequence{DNAAlphabet{2}}("TTAGC")
5nt DNA Sequence:
TTAGC

julia> LongRNASeq("UUANC")
5nt RNA Sequence:
UUANC

julia> LongSequence{RNAAlphabet{2}}("UUAGC")
5nt RNA Sequence:
UUAGC

julia> DNAMer{8}("ATCGATCG")
DNA 8-mer:
ATCGATCG

julia> # The following works, but is not type-stable...

julia> DNAMer("ATCGATCG")
DNA 8-mer:
ATCGATCG

julia> RNAMer{8}("AUCGAUCG")
RNA 8-mer:
AUCGAUCG

julia> # The following works, but is not type-stable...

julia> RNAMer(LongRNASeq("AUCGAUCG"))
RNA 8-mer:
AUCGAUCG

julia> seq = ReferenceSequence("NNCGTATTTTCN")
12nt Reference Sequence:
NNCGTATTTTCN

```

!!! note 
    From version 2.0 onwards, the `convert` methods for converting a
    string or vector of symbols into a sequence type have been removed. These
    convert methods did nothing but pass their arguments to the appropriate
    constructor.
    
    These specific `convert` methods have been removed due to the semantics of
    `convert`: Even though `convert(LongDNASeq, "ATCG")` was previously the same
    as `LongDNASeq("ATCG")`, unlike constructors, `convert` is sometimes
    implicitly called. So it's methods should be restricted to cases that are
    considered safe or unsurprising. `convert` should convert between types that
    represent the same basic kind of thing, like different representations of
    numbers. It is also usually lossless. Not all strings are valid sequences,
    and depending on the sequence type, not all vectors of BioSymbols are valid
    sequences either. A string only represents the "Same kind of thing" as a
    biological sequence in some cases, so implicitly `convert` ing things to a
    sequence type was never safe or unsurprising. These `convert` methods have
    been changed to `Base.parse` methods.

### Constructing sequences from arrays of BioSymbols

Sequences can be constructed using vectors or arrays of a `BioSymbol` type:

```jldoctest
julia> LongDNASeq([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

julia> LongSequence{DNAAlphabet{2}}([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
5nt DNA Sequence:
TTAGC

julia> DNAMer{5}([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
DNA 5-mer:
TTAGC

julia> RNAMer{5}([RNA_U, RNA_U, RNA_A, RNA_G, RNA_C])
RNA 5-mer:
UUAGC

julia> # Works, but is not type-stable

julia> DNAMer([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
DNA 5-mer:
TTAGC

julia> RNAMer([RNA_U, RNA_U, RNA_A, RNA_G, RNA_C])
RNA 5-mer:
UUAGC
```

### Constructing sequences from other sequences

You can create sequences, by concatenating other sequences together:
```jldoctest
julia> LongDNASeq(LongDNASeq("ACGT"), LongDNASeq("NNNN"), LongDNASeq("TGCA"))
12nt DNA Sequence:
ACGTNNNNTGCA

julia> LongDNASeq("ACGT") * LongDNASeq("TGCA")
8nt DNA Sequence:
ACGTTGCA

julia> repeat(LongDNASeq("TA"), 10)
20nt DNA Sequence:
TATATATATATATATATATA

julia> LongDNASeq("TA") ^ 10
20nt DNA Sequence:
TATATATATATATATATATA

```

You can also construct long sequences from kmer sequences, and vice versa:

```jldoctest
julia> m = DNAMer{5}([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
DNA 5-mer:
TTAGC

julia> LongSequence(m)
5nt DNA Sequence:
TTAGC

julia> # round trip from mer to long sequence back to mer.

julia> DNAMer(LongSequence(m))
DNA 5-mer:
TTAGC
```

## Conversion of sequence types

Sometimes you can convert between sequence types without construction / having
to copy data. for example, despite being separate types, `LongDNASeq` and
`LongRNASeq` can freely be converted between efficiently, without copying the
underlying data:
```jldoctest
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> rna = convert(LongRNASeq, dna)
12nt RNA Sequence:
UUANGUAGACCG

julia> dna.data === rna.data  # underlying data are same
true

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
1-element Array{String,1}:
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

julia> char"αβγδϵ"
5char Char Sequence:
αβγδϵ

```

However, it should be noted that by default these sequence literals
allocate the `LongSequence` object before the code containing the sequence
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
a new DNA sequence variable `CTT` is created, and the `A` nucleotide is pushed
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
the sequence will be allocated whilst the code is running, and not before.
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
```@meta
DocTestSetup = quote
    using BioSequences
end
```

So the take home message of sequence literals is this:

Be careful when you are using sequence literals inside of functions, and inside
the bodies of things like for loops. And if you use them and are unsure, use the
 's' and 'd' flags to ensure the behaviour you get is the behaviour you intend.

### Kmer literals

You can create literals for `Mer`s and `BigMer`s as well:

```jldoctest
julia> mer"ATCG"
DNA 4-mer:
ATCG

julia> mer"ATCG"dna
DNA 4-mer:
ATCG

julia> mer"AUCG"rna
RNA 4-mer:
AUCG

```

By using a flag at the end of the literal, you can set whether the kmer should
be a DNA kmer or an RNA kmer. If you don't set the flag, then by default it will
try to make a dna kmer from the string.

Literals for `BigMer`s are also available:

```jldoctest
julia> bigmer"ATCG"
DNA 4-mer:
ATCG

julia> bigmer"ATCG"dna
DNA 4-mer:
ATCG

julia> bigmer"AUCG"rna
RNA 4-mer:
AUCG

```
