using Documenter, BioSequences

DocMeta.setdocmeta!(BioSequences, :DocTestSetup, :(using BioSequences); recursive=true)

makedocs(
    format = Documenter.HTML(),
    sitename = "BioSequences.jl",
    pages = [
        "Home"                           => "index.md",
        "Biological Symbols"             => "symbols.md",
        "BioSequences Types"             => "types.md",
        "Constructing sequences"         => "construction.md",
        "Indexing & modifying sequences" => "transforms.md",
        "Predicates"                     => "predicates.md",
        "Random sequences"               => "random.md",
        "Pattern matching and searching" => "sequence_search.md",
        "Counting"                       => "counting.md",
        "I/O"                            => "io.md",
        "Implementing custom types"      => "interfaces.md",
        "Recipes"                        => "recipes.md",
    ],
    authors = "Sabrina Jaye Ward, Jakob Nissen, D.C.Jones, Kenta Sato, The BioJulia Organisation and other contributors.",
    checkdocs = :all,
)

deploydocs(
    repo = "github.com/BioJulia/BioSequences.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)
