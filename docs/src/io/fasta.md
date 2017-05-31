```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# IO - FASTA formatted files

FASTA is a text-based file format for representing biological sequences.
A FASTA file stores a list of sequence records with name, description, and
sequence.

The template of a sequence record is:

```
>{name} {description}?
{sequence}
```

Here is an example of a chromosomal sequence:

```
>chrI chromosome 1
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC
CACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG
```

## Readers and Writers
The reader and writer for FASTA formatted files, are found within the
`BioSequences.FASTA` module.

```@docs
FASTA.Reader
FASTA.Writer
```

They can be created with IOStreams:

```jlcon
r = FASTA.Reader(open("MyInput.fasta", "r"))
w = FASTA.Writer(open("MyFile.fasta", "w"))
```

Usually sequence records will be read sequentially from a file by iteration.

```jlcon
using BioSequences
reader = FASTA.Reader(open("hg38.fa", "r"))
for record in reader
    # Do something
end
close(reader)
```

But if the FASTA file has an auxiliary index file formatted in fai, the reader
supports random access to FASTA records, which would be useful when accessing
specific parts of a huge genome sequence:

```jlcon
reader = open(FASTAReader, "sacCer.fa", index="sacCer.fa.fai")
chrIV = reader["chrIV"]  # directly read sequences called chrIV.
```

Reading in a sequence from a FASTA formatted file will give you a variable of
type `FASTA.Record`.

```@docs
FASTA.Record
```

Various getters and setters are available for `FASTA.Record`s:

```@docs
FASTA.hasidentifier
FASTA.identifier
FASTA.hasdescription
FASTA.description
FASTA.hassequence
FASTA.sequence(record::FASTA.Record, [part::UnitRange{Int}])
```

To write a `BioSequence` to FASTA file, you first have to create a `FASTA.Record`:

```jlcon
using BioSequences
x = dna"aaaaatttttcccccggggg"
rec = FASTA.Record("MySeq", x)
w = FASTA.Writer(open("MyFile.fasta", "w"))
write(w, rec)
```

As always with julia IO types, remember to close your file readers and writer
after you are finished.
