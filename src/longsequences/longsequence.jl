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
#   data : ....|xxxxx...........|xxXxxxxxxxxxxxxx|............xxxx|....
# offset :                          |<-offset(j)-|
#  width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
#
#  * '.' : unused (4 bits/char)
#  * 'x' : used
#  * 'X' : used and pointed by index `i`

"""
Biological sequence data structure indexed by an alphabet type `A`.
"""
mutable struct LongSequence{A <: Alphabet} <: BioSequence{A}
    data::Vector{UInt64}  # encoded character sequence data
    len::Int

    function LongSequence{A}(data::Vector{UInt64}, len::Int) where {A <: Alphabet}
        return new(data, len)
    end
end

const LongNucleotideSequence = LongSequence{<:NucleicAcidAlphabet}
const LongDNASeq       = LongSequence{DNAAlphabet{4}}
const LongRNASeq       = LongSequence{RNAAlphabet{4}}
const LongAminoAcidSeq = LongSequence{AminoAcidAlphabet}
const LongCharSeq      = LongSequence{CharAlphabet}

###
### Required type traits and methods
###

"Gets the alphabet encoding of a given BioSequence."
BioSymbols.alphabet(::Type{LongSequence{A}}) where {A} = alphabet(A)
Alphabet(::Type{LongSequence{A}}) where {A <: Alphabet} = A()
Base.length(seq::LongSequence) = seq.len
bindata(seq::LongSequence) = seq.data
Base.eltype(::Type{LongSequence{A}}) where {A} = eltype(A)

@inline seq_data_len(s::LongSequence{A}) where A = seq_data_len(A, length(s))
@inline function seq_data_len(::Type{A}, len::Integer) where A <: Alphabet
    return cld(len, div(64, bits_per_symbol(A())))
end

@inline function encoded_data(seq::LongSequence)
    return seq.data
end

include("indexing.jl")
include("constructors.jl")
include("printing.jl")
include("copying.jl")
include("conversion.jl")
include("stringliterals.jl")
include("transformations.jl")
include("operators.jl")
include("counting.jl")
