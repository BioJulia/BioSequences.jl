using Documenter, BioSequences

makedocs(
    format = Documenter.HTML(),
    sitename = "BioSequences.jl",
    pages = [
        "Home"                           => "index.md",
        "Biological Symbols"             => "symbols.md",
        "BioSequences Types"             => "types.md",
        "Constructing sequences"         => "construction.md",
        "Indexing & modifying sequences" => "transforms.md",
        "Random sequences"               => "random.md",
        "Pattern matching and searching" => "search.md",
        "Sequence Demultipllexing"       => "demultiplexer.md"

        #"User Manual" => [
        #    "Sequence types" => [
        #        "Overview" => "user_manual/sequences/biosequence.md",
        #        "Long Sequence" => "user_manual/sequences/generalseq.md",
        #        "Reference Sequences" => "user_manual/sequences/refseq.md",
        #        "Skipmers and kmers" => "user_manual/sequences/skipmers.md"
        #    ],
        #    "Indexing sequences" => "user_manual/indexing.md",
        #    "Searching" => "user_manual/search.md",
        #    "Sequence Composition" => "user_manual/composition.md",
        #    "Demultiplexing" => "user_manual/demultiplexer.md",
        #],
        #"Developer Notes" => [
        #    "Biological symbols" => "dev_docs/symbols.md",
        #    "Biological sequence types" => [
        #        "BioSequence" => "dev_docs/sequences/biosequence.md"
        #    ],
        #    "BitIndex" => "dev_docs/bitindex.md"
        #]
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
