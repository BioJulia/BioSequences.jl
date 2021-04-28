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

struct BitIndex{N,W}
    val::Int64
end

BitsPerSymbol(::BitIndex{N, W}) where {N,W} = BitsPerSymbol{N}()
bits_per_symbol(::BitIndex{N, W}) where {N,W} = N

@inline function bitindex(::BitsPerSymbol{N}, ::Type{W}, i) where {N,W}
    return BitIndex{N, W}((i - 1) << trailing_zeros(N))
end

@inline bitwidth(::Type{W}) where {W<:Unsigned} = 8 * sizeof(W)
@inline bitwidth(::BitIndex{N,W}) where {N,W} = bitwidth(W)

@inline index_shift(i::BitIndex{N,W}) where {N,W} = trailing_zeros(bitwidth(W))
@inline offset_mask(i::BitIndex{N,W}) where {N,W} = UInt8(bitwidth(W)) - 0x01

@inline index(i::BitIndex) = (i.val >> index_shift(i)) + 1
@inline offset(i::BitIndex) = i.val & offset_mask(i)

Base.:+(i::BitIndex{N,W}, n::Integer) where {N,W} = BitIndex{N,W}(i.val + n)
Base.:-(i::BitIndex{N,W}, n::Integer) where {N,W} = BitIndex{N,W}(i.val - n)
Base.:-(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
Base.:(==)(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)
Base.cmp(i1::BitIndex, i2::BitIndex) = cmp(i1.val, i2.val)

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
@inline function extract_encoded_element(bidx::BitIndex{N,W}, data::AbstractArray{W}) where {N,W}
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::BitIndex{N,W}, data::NTuple{n,W}) where {N,n,W}
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> (bitwidth(bidx) - N - offset(bidx))
    return offchunk & bitmask(bidx)
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
bitmask(bidx::BitIndex{N,W}) where {N, W} = bitmask(W, N)
bitmask(n::Integer) = bitmask(UInt64, n)
bitmask(::Type{T}, ::Val{N}) where {T, N} = (one(T) << N) - one(T)


# TODO: Work out places this is used and see if it is really nessecery given the
# bitmask methods above.
# TODO: Resolve this use of bits_per_symbol and A().
bitmask(::A) where {A<:Alphabet} = bitmask(bits_per_symbol(A()))
