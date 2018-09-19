# ShortSequence
# =============
#
# A compact short sequence type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Compact short sequence types
    
Sometimes applications do not use long sequences, but short ones.
Such applications include tools that analyse and assemble kmers generated from
reads, or tools that analyse codons: short 3bp long sequences that correspond to
amino acids when a gene is translated into its protein product.
It is conveinient to represent such sequences using single 'words'.
    
Representation
--------------

Four kinds of nucleotides are encoded as follows:

| nucleotide | binary |
|:----------:|:------:|
|     A      |   00   |
|     C      |   01   |
|     G      |   10   |
|   T / U    |   11   |
    
Nucleic acids are filled from MSBs to LSBs and right-aligned so that all
short sequences are lexicographically ordered.
For example, the memory layout of "TACG" is:

```
64-bit: 0b 00 00 … 00 11 00 01 10
 4-mer:                T  A  C  G
```

With this representation, the unused portion of the unsigned integer, should
be kept blank (all zeros). 
"""
abstract type ShortSequence{N} <: BioSequence end

encoded_data_eltype(::Type{T}) where T <: ShortSequence{128} = UInt128
encoded_data_eltype(::Type{T}) where T <: ShortSequence{64}  = UInt64
encoded_data_eltype(::Type{T}) where T <: ShortSequence{32}  = UInt32
encoded_data_eltype(::Type{T}) where T <: ShortSequence{16}  = UInt16
encoded_data_eltype(::Type{T}) where T <: ShortSequence{8}   = UInt8
encoded_data_eltype(x::ShortSequence) = encoded_data_eltype(typeof(x))

encoded_data(x::ShortSequence) = reinterpret(encoded_data_eltype(x), x)

@inline capacity(::Type{<:ShortSequence{N}}) where N = div(N, 2)
@inline capacity(x::ShortSequence) = capacity(typeof(x))
@inline n_unused(x::ShortSequence) where N = capacity(x) - length(x)

include("indexing.jl")
include("predicates.jl")
include("operations.jl")
include("transformations.jl")

Base.:-(x::ShortSequence, y::Integer) where {T,K} = typeof(x)(encoded_data(x) - y % encoded_data_eltype(x))
Base.:+(x::ShortSequence, y::Integer) where {T,K} = typeof(x)(encoded_data(x) + y % encoded_data_eltype(x))
Base.:+(x::Integer, y::ShortSequence) where {T,K} = y + x

Base.hash(x::ShortSequence, h::UInt) = hash(encoded_data(x), h)
