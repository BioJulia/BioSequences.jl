```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# Biological symbols

The `BioSequences` module reexports the biological symbol (character) types
that are provided by [BioSymbols.jl:](https://biojulia.net/BioSymbols.jl/)

| Type            | Meaning        |
| :-------------- | :------------- |
| `DNA`           | DNA nucleotide |
| `RNA`           | RNA nucleotide |
| `AminoAcid`     | Amino acid     |

These symbols are elements of biological sequences, just as characters are
elements of strings. See sections beginning from
[Sequence Types](@ref) section for details.


## DNA and RNA nucleotides

Set of nucleotide symbols in BioSequences covers IUPAC nucleotide base plus
a gap symbol:

| Symbol | Constant              | Meaning                    |
| :----- | :-------------------- | :------------------------- |
| 'A'    | `DNA_A` / `RNA_A`     | A; Adenine                 |
| 'C'    | `DNA_C` / `RNA_C`     | C; Cytosine                |
| 'G'    | `DNA_G` / `RNA_G`     | G; Guanine                 |
| 'T'    | `DNA_T`               | T; Thymine (DNA only)      |
| 'U'    | `RNA_U`               | U; Uracil (RNA only)       |
| 'M'    | `DNA_M` / `RNA_M`     | A or C                     |
| 'R'    | `DNA_R` / `RNA_R`     | A or G                     |
| 'W'    | `DNA_W` / `RNA_W`     | A or T/U                   |
| 'S'    | `DNA_S` / `RNA_S`     | C or G                     |
| 'Y'    | `DNA_Y` / `RNA_Y`     | C or T/U                   |
| 'K'    | `DNA_K` / `RNA_K`     | G or T/U                   |
| 'V'    | `DNA_V` / `RNA_V`     | A or C or G; not T/U       |
| 'H'    | `DNA_H` / `RNA_H`     | A or C or T; not G         |
| 'D'    | `DNA_D` / `RNA_D`     | A or G or T/U; not C       |
| 'B'    | `DNA_B` / `RNA_B`     | C or G or T/U; not A       |
| 'N'    | `DNA_N` / `RNA_N`     | A or C or G or T/U         |
| '-'    | `DNA_Gap` / `RNA_Gap` | Gap (none of the above)    |

<https://www.bioinformatics.org/sms/iupac.html>

Symbols are accessible as constants with `DNA_` or `RNA_` prefix:
```jldoctest
julia> DNA_A
DNA_A

julia> DNA_T
DNA_T

julia> RNA_U
RNA_U

julia> DNA_Gap
DNA_Gap

julia> typeof(DNA_A)
DNA

julia> typeof(RNA_A)
RNA

```

Symbols can be constructed by converting regular characters:
```jldoctest
julia> convert(DNA, 'C')
DNA_C

julia> convert(DNA, 'C') === DNA_C
true

```


## Amino acids

Set of amino acid symbols also covers IUPAC amino acid symbols plus a gap symbol:

| Symbol       | Constant        | Meaning                     |
| :----------- | :-------------- | :-------------------------- |
| 'A'          | `AA_A`          | Alanine                     |
| 'R'          | `AA_R`          | Arginine                    |
| 'N'          | `AA_N`          | Asparagine                  |
| 'D'          | `AA_D`          | Aspartic acid (Aspartate)   |
| 'C'          | `AA_C`          | Cysteine                    |
| 'Q'          | `AA_Q`          | Glutamine                   |
| 'E'          | `AA_E`          | Glutamic acid (Glutamate)   |
| 'G'          | `AA_G`          | Glycine                     |
| 'H'          | `AA_H`          | Histidine                   |
| 'I'          | `AA_I`          | Isoleucine                  |
| 'L'          | `AA_L`          | Leucine                     |
| 'K'          | `AA_K`          | Lysine                      |
| 'M'          | `AA_M`          | Methionine                  |
| 'F'          | `AA_F`          | Phenylalanine               |
| 'P'          | `AA_P`          | Proline                     |
| 'S'          | `AA_S`          | Serine                      |
| 'T'          | `AA_T`          | Threonine                   |
| 'W'          | `AA_W`          | Tryptophan                  |
| 'Y'          | `AA_Y`          | Tyrosine                    |
| 'V'          | `AA_V`          | Valine                      |
| 'O'          | `AA_O`          | Pyrrolysine                 |
| 'U'          | `AA_U`          | Selenocysteine              |
| 'B'          | `AA_B`          | Aspartic acid or Asparagine |
| 'J'          | `AA_J`          | Leucine or Isoleucine       |
| 'Z'          | `AA_Z`          | Glutamine or Glutamic acid  |
| 'X'          | `AA_X`          | Any amino acid              |
| '*'          | `AA_Term`       | Termination codon           |
| '-'          | `AA_Gap`        | Gap (none of the above)     |

<https://www.bioinformatics.org/sms/iupac.html>

Symbols are accessible as constants with `AA_` prefix:
```jldoctest
julia> AA_A
AA_A

julia> AA_Q
AA_Q

julia> AA_Term
AA_Term

julia> typeof(AA_A)
AminoAcid

```

Symbols can be constructed by converting regular characters:
```jldoctest
julia> convert(AminoAcid, 'A')
AA_A

julia> convert(AminoAcid, 'P') === AA_P
true

```


## Other functions

```@docs
alphabet
gap
iscompatible
isambiguous
```
