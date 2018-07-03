```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# Biological symbols: Developer Information

## Bit encoding of biological symbols

### Nucleic acids (`DNA` and `RNA`)

Every nucleotide is encoded using the lower 4 bits of a byte.
An unambiguous nucleotide has only one set bit and the other bits are unset.
The table below summarises all unambiguous nucleotides and their corresponding bits.
An ambiguous nucleotide is the bitwise OR of unambiguous nucleotides that the
ambiguous nucleotide can take.
For example, `DNA_R` (meaning the nucleotide is either `DNA_A` or `DNA_G`) is
encoded as `0101` because `0101` is the bitwise OR of `0001` (`DNA_A`) and
`0100` (`DNA_G`).
The gap symbol is always `0000`.

|   NucleicAcid    |  Bits  |
|:---------------- |:------ |
| `DNA_A`, `RNA_A` | `0001` |
| `DNA_C`, `RNA_C` | `0010` |
| `DNA_G`, `RNA_G` | `0100` |
| `DNA_T`, `RNA_U` | `1000` |

The next few examples demonstrate some of the bit operations of DNA to illustrate:
```jldoctest
julia> bits(reinterpret(UInt8, DNA_A))
"00000001"

julia> bits(reinterpret(UInt8, DNA_G))
"00000100"

julia> bits(reinterpret(UInt8, DNA_R))
"00000101"

julia> bits(reinterpret(UInt8, DNA_B))
"00001110"

julia> ~DNA_A
DNA_B

julia> DNA_A | DNA_G
DNA_R

julia> DNA_R & DNA_B
DNA_G

```
