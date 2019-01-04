# GeneralSequence
# ==================
#
# A general purpose biological sequence representation.
#
# It is cheap to create a subsequence from a sequence because sequences can
# share underlying data: creating a subsequence is just allocating a new
# sequence object defining a part of the original sequence without copy.
# Destructive operations create a copy of the underlying data if and only if it
# is shared between (sub)sequences. This is often called as copy-on-write
# strategy in computer science and should be transparent to the user.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# About Internals
# ---------------
#
# The `data` field of a `GeneralSequence{A}` object contains binary representation
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
mutable struct GeneralSequence{A <: Alphabet} <: BioSequence{A}
    data::Vector{UInt64}  # encoded character sequence data
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
    shared::Bool          # true if and only if `data` is shared between sequences

    function GeneralSequence{A}(data::Vector{UInt64},
                            part::UnitRange{Int},
                            shared::Bool) where A
        return new(data, part, shared)
    end
end

const DNASequence       = GeneralSequence{DNAAlphabet{4}}
const RNASequence       = GeneralSequence{RNAAlphabet{4}}
const AminoAcidSequence = GeneralSequence{AminoAcidAlphabet}
const CharSequence      = GeneralSequence{CharAlphabet}


# Required type traits and methods
# ================================

"Gets the alphabet encoding of a given BioSequence."
BioSymbols.alphabet(::Type{GeneralSequence{A}}) where {A} = alphabet(A)
Alphabet(::Type{GeneralSequence{A}}) where {A <: Alphabet} = A()
Base.length(seq::GeneralSequence) = length(seq.part)
bindata(seq::GeneralSequence) = seq.data
Base.eltype(::Type{GeneralSequence{A}}) where {A} = eltype(A)

@inbounds function seq_data_len(::Type{A}, len::Integer) where A <: Alphabet
    return cld(len, div(64, bits_per_symbol(A())))
end

@inbounds function encoded_data(seq::GeneralSequence)
    return seq.data
end


# Replace a GeneralSequence's data with a copy, copying only what's needed.
# The user should never need to call this, as it has no outward effect on the
# sequence.
function orphan!(seq::GeneralSequence{A},
		 size::Integer = length(seq),
		 force::Bool = false) where {A}

    if !seq.shared && !force
        return seq
    end

    j, r = bitindex(seq, 1)
    data = Vector{UInt64}(undef, seq_data_len(A, size))

    if !isempty(seq) && !isempty(data)
        x = seq.data[j] >> r
        m = index(bitindex(seq, lastindex(seq))) - j + 1
        l = min(lastindex(data), m)
        @inbounds @simd for i in 1:l-1
            y = seq.data[j + i]
            data[i] = x | y << (64 - r)
            x = y >> r
        end
        if m <= l
            data[l] = x
        else
            y = seq.data[j + l]
            data[l] = x | y << (64 - r)
        end
    end

    seq.data = data
    seq.part = 1:length(seq)
    seq.shared = false
    return seq
end

include("indexing.jl")
include("constructors.jl")
include("iteration.jl")
include("copying.jl")
include("conversion.jl")
include("stringliterals.jl")
include("transformations.jl")
include("operators.jl")
