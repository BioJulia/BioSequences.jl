```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# IO - 2bit formatted files

2bit is a binary file format designed for storing a genome consists of multiple
chromosomal sequences.
The reading speed is often an order of magnitude faster than that of FASTA and
the file size is smaller.
However, since the .2bit file format is specialized for genomic sequences, it
cannot store either RNA or amino acid sequences.

## Readers and Writers

The reader and writer for 2bit formatted files, are found within the
`BioSequences.TwoBit` module.

```@docs
TwoBit.Reader
TwoBit.Writer
```

The 2bit reader supports random access using an index included in the header
section of a .2bit file:

```jlcon
reader = TwoBit.Reader(open("sacCer.2bit", "r"))
chrIV = reader["chrIV"] # directly read chromosome 4
```

If you want to know the names of the sequences available in the file,
you can use the `seqnames` method on the reader.

```jlcon
seqnames(reader)
```

Reading from a `TwoBit.Reader` will yield a `TwoBit.Record` type variable:

```@docs
TwoBit.Record
```

To write a sequence to a TwoBit file, first a record must be created.

```@docs
TwoBit.Record(name::AbstractString, seq::BioSequences.Sequence, masks::Union{Vector{UnitRange{Int}}, Nothing})
```
