using Documenter, BioSequences

makedocs(
    format = Documenter.HTML(),
    sitename = "BioSequences.jl",
    pages = [
        "Home" => "index.md",
        "Biological Symbols" => "symbols.md",
        "Sequence types" => [
            "Overview" => "sequences/biosequences.md",
            "BioSequence" => "sequences/bioseq.md",
            "Reference Sequences" => "sequences/refseq.md",
            "Nucleic acid k-mers" => "sequences/kmer.md"
        ],
        "IO" => [
            "FASTA formatted files" => "io/fasta.md",
            "FASTQ formatted files" => "io/fastq.md",
            "2bit formatted files" => "io/twobit.md"
        ],
        "Searching" => "searching.md",
        "Sequence Composition" => "composition.md",
        "Demultiplexing" => "demultiplexer.md",
        "Contributing" => "contributing.md"
    ],
    authors = "Ben J. Ward, D.C.Jones, Kenta Sato, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/BioSequences.jl.git",
    julia = "1.0",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
