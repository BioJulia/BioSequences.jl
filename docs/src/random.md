# Generating random sequences

## Long sequences

You can generate random long sequences using the `randdna` function and the
`Sampler`'s implemented in BioSequences:

```@docs
randdna
randdnaseq
randrnaseq
randaaseq
SamplerUniform
SamplerWeighted
```

## Kmer sequences

You can make random `Mer` quite simply using `Base.rand`:

```jldoctest
julia> rand(DNAMer{7})
DNA 7-mer:
GAATTGC

julia> rand(RNAMer{8})
RNA 8-mer:
UUUCAGUA

julia> rand(BigDNAMer{63})
DNA 63-mer:
TAGATTATCCAGGGCCTTTGACTGGCTCGTATCGACGTAAGTGCCCCGCGGGACATAGGCTGC

```