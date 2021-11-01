```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# Generating random sequences

## Long sequences

You can generate random long sequences using the `randdna` function and the
`Sampler`'s implemented in BioSequences:

```@docs
randseq
randdnaseq
randrnaseq
randaaseq
SamplerUniform
SamplerWeighted
```

## Kmer sequences

You can make random `Mer` quite simply using `Base.rand`:

```@repl
using BioSequences # hide
rand(DNAMer{7})
rand(RNAMer{8})
rand(BigDNAMer{63})
```
