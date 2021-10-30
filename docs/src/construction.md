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

You can convert between sequence types, if the sequences are compatible - that is, if the source sequence does not contain symbols that are un-encodable by the destination type.
```jldoctest
julia> dna = dna"TTACGTAGACCG"
12nt DNA Sequence:
TTACGTAGACCG

julia> dna2 = convert(LongDNA{2}, dna)
12nt DNA Sequence:
TTACGTAGACCG
```

DNA/RNA are special in that they can be converted to each other, despite containing distinct symbols.
When doing so, `DNA_T` is converted to `RNA_U` and vice versa.
```jldoctest
julia> convert(LongRNA{2}, dna"TAGCTAGG")
8nt RNA Sequence:
UAGCUAGG
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

## Comparison to other sequence types
Following Base standards, BioSequences do not compare equal to other containers even if they have the same elements.
To e.g. compare a BioSequence with a vector of DNA, compare the elements themselves:
```jldoctest
julia> seq = dna"GAGCTGA"; vec = collect(seq);

julia> seq == vec, isequal(seq, vec)
(false, false)

julia> length(seq) == length(vec) && all(i == j for (i, j) in zip(seq, vec))
true 
```