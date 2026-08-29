```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# Generating random sequences

## Long sequences

You can generate random long sequences using the `randdnaseq` function and the
samplers implemented in BioSequences:

```@docs
randseq
randdnaseq
randrnaseq
randaaseq
SamplerUniform
SamplerWeighted
```
