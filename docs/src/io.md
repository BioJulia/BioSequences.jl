# I/O for sequencing file formats

Versions of BioSequences prior to v2.0 provided a FASTA, FASTQ, and 2Bit
submodule for working with formatted sequence files.

After version v2.0, in order to neatly separate concerns, these submodules were
removed.

Instead there will now be dedicated BioJulia packages for each format. Each
of these will be compatible with BioSequences.

A list of all of the different formats and packages is provided below to help
you find them quickly.

| Format | Package                                            |
|:------ |:------------------------------------------------   |
| FASTA  | [FASTX.jl](https://github.com/BioJulia/FASTX.jl)   |
| FASTQ  | [FASTX.jl](https://github.com/BioJulia/FASTX.jl)   |
| 2Bit   | [TwoBit.jl](https://github.com/BioJulia/TwoBit.jl) |