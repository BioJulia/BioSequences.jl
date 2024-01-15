```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
    using BioSymbols
end
```

## Recipes
This page provides tested example code to solve various common problems using
BioSequences.

### One-hot encoding biosequences
The types `DNA`, `RNA` and `AminoAcid` expose a binary representation through
the exported function `BioSymbols.compatbits`, which is a one-hot encoding of:

```jldoctest
julia> using BioSymbols

julia> compatbits(DNA_W)
0x09

julia> compatbits(AA_J)
0x00000600
```

Each set bit in the encoding corresponds to a compatible unambiguous symbol.
For example, for `RNA`, the four lower bits encode A, C, G, and U, in order.
Hence, the symbol `D`, which is short for A, G or U, is encoded as
`0x01 | 0x04 | 0x08 == 0x0d`:

```jldoctest
julia> compatbits(RNA_D)
0x0d

julia> compatbits(RNA_A) | compatbits(DNA_G) | compatbits(RNA_U)
0x0d
```

Using this, we can construct a function to one-hot encode sequences - in this
example, nucleic acid sequences:
```jldoctest
function one_hot(s::NucSeq)
    M = falses(4, length(s))
    for (i, s) in enumerate(s)
        bits = compatbits(s)
        while !iszero(bits)
            M[trailing_zeros(bits) + 1, i] = true
            bits &= bits - one(bits) # clear lowest bit
        end
    end
    M
end

one_hot(dna"TGNTKCTW-T")

# output

4×10 BitMatrix:
 0  0  1  0  0  0  0  1  0  0
 0  0  1  0  0  1  0  0  0  0
 0  1  1  0  1  0  0  0  0  0
 1  0  1  1  1  0  1  1  0  1
```


