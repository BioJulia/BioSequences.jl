# BitIndex
# --------
#
# Utils for indexing bits in a vector of unsigned integers (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md
#
#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

"""
    BitIndex

`BitIndex` is an internal type used in BioSequences.It contains
a bit offset. For biosequences with an internal array of coding units,
it can be used to obtain the array index and element bit offset.

Useful methods:
* bitindex(::BioSequence, ::Int)
* index(::BitIndex)
* offset(::BitIndex)
* nextposition / prevposition(::BitIndex)
* extract_encoded_element(::BitIndex, ::Union{Array, Tuple})
"""
struct BitIndex{N, W}
    val::UInt64
end

BitsPerSymbol(::BitIndex{N, W}) where {N,W} = BitsPerSymbol{N}()
bits_per_symbol(::BitIndex{N, W}) where {N,W} = N

@inline function bitindex(::BitsPerSymbol{N}, ::Type{W}, i) where {N,W}
    return BitIndex{N, W}((i - 1) << trailing_zeros(N))
end

bitwidth(::Type{W}) where W = 8*sizeof(W)
@inline index_shift(::BitIndex{N,W}) where {N,W} = trailing_zeros(bitwidth(W))
@inline offset_mask(::BitIndex{N,W}) where {N,W} = UInt8(bitwidth(W)) - 0x01

@inline index(i::BitIndex) = (i.val >> index_shift(i)) + 1
@inline offset(i::BitIndex) = i.val & offset_mask(i)

Base.:+(i::BitIndex{N,W}, n::Integer) where {N,W} = BitIndex{N,W}(i.val + n)
Base.:-(i::BitIndex{N,W}, n::Integer) where {N,W} = BitIndex{N,W}(i.val - n)
Base.:-(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
Base.:(==)(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)

@inline function nextposition(i::BitIndex{N,W}) where {N,W}
    return i + N
end

@inline function prevposition(i::BitIndex{N,W}) where {N,W}
    return i - N
end

function Base.iterate(i::BitIndex, s = 1)
    if s == 1
        return (index(i), 2)
    elseif s == 2
        return (offset(i), 3)
    else
        return nothing
    end
end

Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::BitIndex{N,W}, data::Union{AbstractArray{W}, Tuple{Vararg{W}}}) where {N,W}
    chunk = (@inbounds data[index(bidx)]) >> offset(bidx)
    return chunk & bitmask(bidx)
end

# Create a bit mask that fills least significant `n` bits (`n` must be a
# non-negative integer).
"Create a bit mask covering the least significant `n` bits."
function bitmask(::Type{T}, n::Integer) where {T}
    topshift = 8 * sizeof(T) - 1
    return ifelse(n > topshift, typemax(T), one(T) << (n & topshift) - one(T))
end

# Create a bit mask filling least significant N bits.
# This is used in the extract_encoded_element function.
bitmask(::BitIndex{N,W}) where {N, W} = bitmask(W, N)
bitmask(n::Integer) = bitmask(UInt64, n)
