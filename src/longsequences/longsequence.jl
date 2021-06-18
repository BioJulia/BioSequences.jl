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

`LongSequence` is the default mutable, variable-length `BioSequence`.
It is suitable for biological sequences whose length is either not
known at compile time, or that is larger than about 100 symbols.

See also: [`Mer`](@ref)
"""
mutable struct LongSequence{A <: Alphabet} <: BioSequence{A}
    data::Vector{UInt64}  # encoded character sequence data
    len::UInt

    function LongSequence{A}(data::Vector{UInt64}, len::UInt) where {A <: Alphabet}
        new{A}(data, len)
    end
end

const LongNucleotideSequence = LongSequence{<:NucleicAcidAlphabet}
const LongDNASeq       = LongSequence{DNAAlphabet{4}}
const LongRNASeq       = LongSequence{RNAAlphabet{4}}
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