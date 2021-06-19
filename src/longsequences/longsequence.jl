###
### LongSequence
###
###
### A general purpose biological sequence representation.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# About Internals
# ---------------
#
# The `data` field of a `LongSequence{A}` object contains binary representation
# of a biological character sequence. Each character is encoded with an encoder
# corresponding to the alphabet `A` and compactly packed into `data`. To extract
# a character from a sequence, you should decode this binary sequence with a
# decoder that is a pair of the encoder. The length of encoded binary bits is
# fixed, and hence a character at arbitrary position can be extracted in a
# constant time. To know the exact location of a character at a position, you
# can use the `bitindex(seq, i)` function, which returns a pair of element's index
# containing binary bits and bits' offset. As a whole, character extraction
# `seq[i]` can be written as:
#
#     j = bitindex(seq, i)
#     decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
#
#  index :        index(j) - 1       index(j)       index(j) + 1
#   data :     |xxxxxxxxxxxxxxxx|xxXxxxxxxxxxxxxx|............xxxx|....
# offset :                          |<-offset(j)-|
#  width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
#
#  * '.' : unused (4 bits/char)
#  * 'x' : used
#  * 'X' : used and pointed by index `i`

"""
    LongSequence{A <: Alphabet}

Many genomics scripts and tools benefit from an efficient general purpose
sequence type that allows you to create and edit sequences.

`LongSequence` is the default mutable, variable-length `BioSequence`.
It is suitable for biological sequences whose length is either not
known at compile time, or that is larger than about 100 symbols.

`LongSequence{A<:Alphabet} <: BioSequence{A}` is parameterized by a concrete
`Alphabet` type `A` that defines the domain (or set) of biological symbols
permitted.

As the [`BioSequence`](@ref) interface definition implies, `LongSequence`s store
the biological symbol elements that they contain in a succinct encoded form that
permits many operations to be done in an efficient bit-parallel manner. As per
the interface of [`BioSequence`](@ref), the [`Alphabet`](@ref) determines how
an element is encoded or decoded when it is inserted or extracted from the
sequence.

For example, [`AminoAcidAlphabet`](@ref) is associated with `AminoAcid` and hence
an object of the `LongSequence{AminoAcidAlphabet}` type represents a sequence of
amino acids.

Symbols from multiple alphabets can't be intermixed in one sequence type.

The following table summarizes common LongSequence types that have been given
aliases for convenience.

| Type                                | Symbol type | Type alias         |
| :---------------------------------- | :---------- | :----------------- |
| `LongSequence{DNAAlphabet{4}}`      | `DNA`       | `LongDNASeq`       |
| `LongSequence{RNAAlphabet{4}}`      | `RNA`       | `LongRNASeq`       |
| `LongSequence{AminoAcidAlphabet}`   | `AminoAcid` | `LongAminoAcidSeq` |

The `LongDNASeq` and `LongRNASeq` aliases use a DNAAlphabet{4}.

`DNAAlphabet{4}` permits ambiguous nucleotides, and a sequence must use at least
4 bits to internally store each element (and indeed `LongSequence` does).

If you are sure that you are working with sequences with no ambiguous nucleotides,
you can use `LongSeqeunces` parameterised with `DNAAlphabet{2}` instead.

`DNAAlphabet{2}` is an alphabet that uses two bits per base and limits to only
unambiguous nucleotide symbols (A,C,G,T).

Changing this single parameter, is all you need to do in order to benefit from memory savings.
Some computations that use bitwise operations will also be dramatically faster.

The same applies with `LongSeqeunce{RNAAlphabet{4}}`, simply replace the alphabet
parameter with `RNAAlphabet{2}` in order to benefit.
"""
mutable struct LongSequence{A <: Alphabet} <: BioSequence{A}
    data::Vector{UInt64}  # encoded character sequence data
    len::UInt

    function LongSequence{A}(data::Vector{UInt64}, len::UInt) where {A <: Alphabet}
        new{A}(data, len)
    end
end

const LongNucleotideSequence = LongSequence{<:NucleicAcidAlphabet}

"An alias for LongSequence{DNAAlphabet{4}}"
const LongDNASeq       = LongSequence{DNAAlphabet{4}}

"An alias for LongSequence{RNAAlphabet{4}}"
const LongRNASeq       = LongSequence{RNAAlphabet{4}}

"An alias for LongSequence{AminoAcidAlphabet}"
const LongAminoAcidSeq = LongSequence{AminoAcidAlphabet}

# Basic attributes
Base.length(seq::LongSequence) = seq.len % Int
encoded_data_eltype(::Type{<:LongSequence}) = UInt64
Base.copy(x::LongSequence) = typeof(x)(copy(x.data), x.len)

# Derived basic attributes
symbols_per_data_element(x::LongSequence) = div(64, bits_per_symbol(Alphabet(x)))

include("seqview.jl")
include("indexing.jl")
include("constructors.jl")
include("conversion.jl")
include("copying.jl")
include("stringliterals.jl")
include("transformations.jl")
include("operators.jl")
include("counting.jl")